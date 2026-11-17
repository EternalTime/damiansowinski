# Chan, Horne and Mann's charged dilaton black holes

Static, spherically symmetric charged black holes of Einstein-Maxwell-dilaton gravity whose region far away is neither flat nor de Sitter nor anti-de Sitter space, with signature $(-,+,+,+)$.
This note records the sources of the eight charts, why the first chart is written in a length $r_0$ and a pure number $h$, the checks each chart passes before it is printed, and what the drawings draw.
Chan, Horne and Mann's article (Nucl. Phys. B 447, 441, 1995, gr-qc/9502042), Clément and Leygnac's (Phys. Rev. D 70, 084018, 2004, gr-qc/0405034) and Clément, Gal'tsov and Leygnac's (Phys. Rev. D 67, 024012, 2003, hep-th/0208225) were read in full.

## Step 1: the action and the field equations

Chan, Horne and Mann extremize $S = \int d^4x\sqrt{-g}\,(\mathcal{R} - 2(\nabla\phi)^2 - V(\phi) - e^{-2a\phi}F^2)$, their (1.1), in units $G = c = 1$.
The field equations are their (2.1): $R_{ab} = 2\partial_a\phi\,\partial_b\phi + \tfrac{1}{2}g_{ab}V + 2e^{-2a\phi}(F_{ac}F_b{}^c - \tfrac{1}{4}g_{ab}F^2)$, $\partial_a(\sqrt{-g}\,e^{-2a\phi}F^{ab}) = 0$ and $\nabla^2\phi = \tfrac{1}{4}\,dV/d\phi - \tfrac{a}{2}e^{-2a\phi}F^2$.
The electric field of their (2.3) is $F_{tr} = Qe^{2a\phi}/R^2$ in the chart of their (2.2), $ds^2 = -U\,dt^2 + dr^2/U + R^2d\Omega^2$, and the magnetic hole has $F_{\theta\phi} = Q\sin\theta$ with $\phi \to -\phi$, as their section 2 says.
Every chart's check in `dilaton_unusual_asymptotics.py` holds these three equations at points, in forty digits, from the chart's own Ricci tensor, for the electric and the magnetic hole where the chart is of the Einstein metric.
The Ricci tensor is the standard contraction, $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$.

## Step 2: the first chart, and why it is written in $r_0$ and $h$

Their solution (3.7) with no potential is $U = r^{2/(1+a^2)}(1 - r_h/r)/\gamma^2$ and $R = \gamma r^{a^2/(1+a^2)}$, with $r_h = 2(1 + a^2)M/a^2$ their (3.8) for the quasilocal mass $M$ of their (2.8).
The constant $\gamma$ carries the dimension of a length to the power $1/(1 + a^2)$.
The checker reads a power whose exponent holds a parameter with its base in a fixed unit, so $r^{2/(1+a^2)}/\gamma^2$ and $\gamma^2r^{2a^2/(1+a^2)}$ cannot both balance with $\gamma$ a parameter of its own.
Clément and Leygnac write the same metric, their (2.8), with a length $r_0$, $r_0^{1+\gamma_{CL}} = \gamma^2$ for their $\gamma_{CL} = (1 - a^2)/(1 + a^2)$, so $\gamma^2 = r_0^{2/(1+a^2)}$, and the chart is written in it.
The chart then names the pure number $h = (r/r_0)^{2/(1+a^2)}$, with $U = h(1 - r_h/r)$ and $R^2 = r^2/h$.
$h$ is held as a function of $r$ with the slope $\partial_rh = 2h/((1 + a^2)r)$, `vm.HELD` and `vm.RATES`, so every value is a rational function of $r$, $h$, $a$ and $r_h$.
Printed with the power written out, the step that factors each value did not finish in ten minutes on 17 November 2026, since sympy factors poorly around exponents in $1/(1 + a^2)$; held, the chart prints in two seconds.
The check holds the chart to their (3.7) with $\gamma^2 = r_0^{2/(1+a^2)}$, to $F_{tr} = 1/((1 + a^2)Q)$, which their (2.3) and (3.7c) give, and to their Ricci scalar (3.9), at couplings $1/2$, $1$, $7/5$ and $3$.
No value of the metric holds $Q$, which enters only through the dilaton, $e^{2a\phi} = R^2/((1 + a^2)Q^2)$.

## Step 3: the charts at $a = 1$

The areal chart is their (3.2), the first chart at $a = 1$ with the areal radius $\gamma r^{1/2}$ as the coordinate, written with $\gamma^2 = r_0$ and $4M = r_h$; the check holds it to their (3.2), to the first chart carried along $r_1 = r^2/r_0$, and to the field equations with $F_{tr} = r/(r_0Q)$, which the first chart's field carried along the map gives.
The electric string metric is their (3.5), $e^{2\phi}g$ in their rescaled radius.
They print the rescaling as $r \to \sqrt{2}\,Qr/\gamma$, but the metric they print follows from $r \to \sqrt{2}\,Qr/\gamma^2$ alone: with $e^{2\phi} = \gamma^2r/(2Q^2)$ from their (3.1c), $\gamma^2r/(2Q^2)$ times their (3.1) in $r = \sqrt{2}\,Q\rho/\gamma^2$ is their (3.5) in $\rho$, and in $r = \sqrt{2}\,Q\rho/\gamma$ it is not.
The check uses the map that gives their (3.5), and holds the chart to $e^{2\phi}$ times the first chart carried along it.
The magnetic string metric is their (3.6), $e^{-2\phi}$ times the first chart at $a = 1$ in the same $r$, with $e^{-2\phi} = 2Q^2/(r_0r)$; it is Witten's black hole in two dimensions times a sphere of radius $\sqrt{2}\,Q$, as they and Giddings and Strominger's (2.13) say.
The chart of no mass is Clément, Gal'tsov and Leygnac's (2.10) and (2.11), with their signature $(+,-,-,-)$ reversed.
Their $r_0 = \sqrt{2}\,\bar{Q}$ of (2.16) is this $r_0$, and their $\alpha^{-2} = 4\sqrt{2}\,\bar{Q}$ is $4r_0$, so $\xi = 2\sqrt{r_0r}\cosh(t/(2r_0))$ and $\eta = 2\sqrt{r_0r}\sinh(t/(2r_0))$ for the $t$ and $r$ of the first chart at $r_h = 0$, and $dt\wedge dr = d\eta\wedge d\xi$; the check holds the chart to the first carried along that map and to the field equations with $F_{\eta\xi} = 1/(2Q)$.
The plane of $\eta$ and $\xi$ is flat and the spheres have the radius $\sqrt{\xi^2 - \eta^2}/2$.

## Step 4: the three families with a Liouville potential

Their (4.1) is $V = 2\Lambda e^{2b\phi}$, and their (4.2), (4.6) and (4.9) are the families at $b = 1/a$, $-a$ and $-1/a$.
(4.2) has $N = 1/(1 + a^2)$, so its own $\gamma^2 = r_0^{2a^2/(1+a^2)}$ and its pure number is $k = (r/r_0)^{2a^2/(1+a^2)}$, with $U_1 = (1 + a^2)k(1/(1 - a^2) - 2M/r + q^2/r^2)$ and $R^2 = r^2/k$ for $q^2 = Q^2e^{2a\phi_0}/\gamma^2$; its horizons are their (4.3).
(4.6) and (4.9) have the first chart's $N$ and $h$: $U_2 = \sigma h(1 - r_h/r)$ for $\sigma = (1 + a^2)(2Q^2e^{2a\phi_0}/\gamma^2 - 1)/(1 - a^2)$, with $r_h$ their (4.7), and $U_3 = h(1 - r_h/r) + \lambda r^2/h$ for $\lambda = (1 + a^2)^{(2a^2 + 1)/a^2}(Q/r_0)^{2/a^2}\Lambda/((1 - 3a^2)a^2)$, their $r^{-2(1-a^2)/(1+a^2)}$ term written through $h$.
At $\sigma = 1$, which is $\Lambda = 0$ in their (4.6d), and at $\lambda = 0$ both are the first chart, as their text says of each.
The check holds each family to their own form at the constants of each point, with the dilaton's $\phi_0$ or the charge $Q$ chosen there and $\Lambda$ from their (4.2d) and (4.6d) or from $\lambda$, and to the field equations with the potential, at three couplings each away from $a^2 = 1$ and $a^2 = 1/3$.
At $a = 0$ the family $b = 1/a$ is Reissner and Nordström's black hole, as their text says.

## Step 5: the drawings

Every drawing that holds the hole embedded is drawn at $a = 1$ and $r_0 = r_h = 1$, Clément, Gal'tsov and Leygnac's linear dilaton black hole, where $h = r$ and $U = r - 1$.
The tortoise coordinate there is $r_* = \ln|r - 1|$, which grows without bound, so a ray takes an infinite $t$ to reach $r \to \infty$ and the conformal diagram is Schwarzschild's, Kruskal's square with $\kappa = 1/2$, $UV = -(r - 1)$ outside and $UV = 1 - r$ inside.
At $a^2 = 1/3$, the value Chan, Horne and Mann name in their section 3.2, $U = r^{1/2}(r - 1)$ and $r_* = \ln|(\sqrt{r} - 1)/(\sqrt{r} + 1)|$, which tends to $0$: a ray reaches $r \to \infty$ at a finite $t$, and $r \to \infty$ is the hyperbola $UV = -1$, the vertical line $X = \pi/2$ of the conformal diagram, a timelike boundary, with the singularity on $UV = 1$ as before.
The spacetime diagram of the first chart is drawn at both couplings; the Liouville families at $a = 1/2$, $M = 4r_0/3$, $q = r_0$ (horizons at $r_0/2$ and $3r_0/2$), at $a = 1/2$, $\sigma = 2$, and at $a = 2$, $\lambda = -1/(10r_h^2)$, where an event horizon at $1.13\,r_h$ and a cosmological horizon at $5.82\,r_h$ bound a static region; the string metrics at $Q = r_0/\sqrt{2}$, which puts both of their horizons at $r = r_h$.
The embedding diagram draws the equatorial plane of $t = 0$ at $a = 1$: its metric is $dr^2/(r - 1) + r\,d\phi^2$, the circle at $r$ has the radius $\rho = \sqrt{r}$, and $(dz/d\rho)^2 = (3r + 1)/(r - 1)$, which tends to $3$, so far away the slice is a cone, and it has a throat at the bifurcation sphere $r = r_h$, where it runs on into the other exterior.
That moment is marked on the first chart's view at $a = 1$, on the areal chart's view, where it runs from $\sqrt{r_0r_h}$ to $\sqrt{6r_0r_h}$, and on the conformal view at $a = 1$.
`exceptions/dilaton_unusual_asymptotics.json` names the drawings on which it does not lie: the other members of the family, and the string metrics, which measure the manifold another way.
The curvature figure plots the stored $K = (11r^2 + 6r\,r_h + 3r_h^2)/(4r_0^2r^4)$ and $R = (r - r_h)/(2r_0r^2)$ of the first chart at $a = 1$, both falling to zero far away, as $11/(4r_0^2r^2)$ and $1/(2r_0r)$, more slowly than Schwarzschild's $K = 12r_s^2/r^6$.
