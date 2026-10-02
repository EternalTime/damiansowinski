# Einstein, Infeld and Hoffmann's field of many bodies

Why each chart is the one published, what the checker counts as an order when the small quantities are potentials, and what `print_charts.py`, `null_rays.py` and `embedding.py` check before they write.

## Step 1: the charts and their sources

The paper is A. Einstein, L. Infeld and B. Hoffmann, "The Gravitational Equations and the Problem of Motion", Annals of Mathematics 39, 65-100 (1938), read in the facsimile of Blum and Rickles's sourcebook, doi:10.34663/9783945561317-17.
Its signature is $(+,-,-,-)$ and its masses are lengths, $G = c = 1$.
Its coordinate conditions are (1,25) and (1,26), $\gamma_{0s|s} - \gamma_{00|0} = 0$ and $\gamma_{ms|s} = 0$, with $\gamma_{\mu\nu}$ the trace-reversed $h_{\mu\nu}$; the first is the time component of the harmonic condition, and the second differs from its space components by $\gamma_{m0|0}$, which is of fourth order, so to the order of the post-Newtonian metric the coordinates are harmonic.
Their footnote 3 says the Lorentz invariant condition "might seem rather more natural" and that theirs makes the calculation simpler.

The field is their (9.8), $h_{00} = \varphi$, $h_{0n} = \gamma_{0n} = \sum 4m\dot\xi^n/r$ and $h_{mn} = \delta_{mn}\varphi$ with $\varphi = \sum(-2m/r)$, and at the next order $h_{00} = \tfrac{1}{2}\gamma_{00} + \tfrac{1}{2}\gamma_{ll}$ from (12.1) with (12.3) and (12.4), and (13.10) with (14.3).
For two bodies that is $h_{00} = \tfrac{1}{2}\varphi^2 - 3m_1v_1^2/r_1 - 3m_2v_2^2/r_2 + (2m_1m_2/d)(1/r_1 + 1/r_2) - m_1\partial_t^2r_1 - m_2\partial_t^2r_2$, with $d$ the distance between them.
In the signature $(-,+,+,+)$, with $U = \sum_a m_a/r_a$, $V_i = \sum_a m_av_a^i/r_a$, $\chi = \sum_a m_ar_a$ and $\psi = \sum_a (m_a/r_a)(3v_a^2/2 - \sum_{b \neq a} m_b/r_{ab})$, it is the harmonic chart:

$g_{00} = -1 + 2U - 2U^2 + 2\psi + \partial_t^2\chi$, $g_{0i} = -4V_i$, $g_{ij} = (1 + 2U)\delta_{ij}$.

Since $\partial_t^2 r_a = (v_a^2 - (\mathbf{n}_a\cdot\mathbf{v}_a)^2)/r_a - \mathbf{n}_a\cdot\mathbf{a}_a$, with the accelerations written by Newton's law this is equation (7.2) of Blanchet, Faye and Ponsot (1998) through its terms in $1/c^4$, $1/c^3$ and $1/c^2$, and Alvi's (2000) near zone metric; `_tools/test_eih_many_bodies.py` holds the published chart to both at a point.

The standard chart is the standard post-Newtonian gauge, Box 2 of Will (2014) at general relativity's values: $g_{00} = -1 + 2U - 2U^2 + 4\Phi_1 + 4\Phi_2$, $g_{0i} = -\tfrac{7}{2}V_i - \tfrac{1}{2}W_i$ and $g_{ij} = (1 + 2U)\delta_{ij}$.
Will's density is the rest mass density in the local frame, $\rho^*(1 - v^2/2 - 3U)$ in terms of the conserved one, so for point masses of constant $m_a$ his $2U + 4\Phi_1 + 4\Phi_2$ is $2U + 2\psi$ here.
For a point mass $W_i - V_i = \partial_t\partial_i(m r)$, so $g_{0i} = -4V_i - \tfrac{1}{2}\partial_t\partial_i\chi$.
The two charts differ by the time alone, $x^0 \to x^0 + \tfrac{1}{2}\partial_t\chi$ from the standard time to the harmonic one, which `eih_check` holds: the harmonic chart pulled back along it is the standard chart, $g_{tt}$ through the fourth order, $g_{ti}$ through the third and $g_{ij}$ through the second.

Both charts are written in the potentials, left free, so they hold for any number of bodies on any paths, and no component assumes a field equation.
The parameters state the equations the potentials solve between the bodies: Laplace's for $U$, $\psi$ and each $V_i$, $\nabla^2\chi = 2U$, and $\partial_tU + \partial_iV_i = 0$.
`eih_between_the_bodies` writes those into a value, and `eih_check` holds every Ricci component kept to vanishing there, and $R_{tt}$, $R_{tx}$ and $R_{xx}$ to not vanishing for free potentials.
The Ricci tensor vanishes for every path of the bodies: the vacuum equations at this order do not fix the motion, which comes from the next order, as the paper's section 11 derives Newton's law from the fourth order of $\Lambda_{mn}$.

A chart of a circular binary with the distances $r_1$ and $r_2$ written out, in corotating coordinates as Alvi has it, was tried first.
With two radicals in every component the checker's geometry did not finish its Riemann tensor in ten minutes, for equal masses or unequal, so that chart is not published; the test holds the harmonic chart to Alvi's metric at a point.

## Step 2: what a post-Newtonian order is when the small quantities are potentials

`ORDERS` gives each potential its order in powers of $v/c$, `({"U": 2, "chi": 2, "V_1": 3, "V_2": 3, "V_3": 3, "psi": 4}, 2, "t")`.
Until this spacetime an order could be declared only for a constant; `Reader` now takes a declared function of the coordinates in a system that names its slow time, and `Reader._scaled` scales the function and each of its derivatives by its power of the one small symbol, a derivative along the time one power more each time it is taken, since $\partial_t \sim (v/c)\,\partial_x$ for sources that move slowly.
That is the paper's own counting, its (3,10) to (3,12): the auxiliary time $\tau = \lambda x^0$, with $\gamma_{00} \sim \lambda^2$, $\gamma_{0n} \sim \lambda^3$ and $\gamma_{mn} \sim \lambda^4$.
Each component is then kept as far as `Reader.kept_order` says, as for `ppn_metric`.
`Orders` in `_tools/test_eih_many_bodies.py` holds the scaling, and holds the rule sound with time dependence: an arbitrary function of all four coordinates in each metric component at the first order left open moves no tensor handed over.

A subscript on a declared function has to be a numeral for `\partial` to bind to it, so the components of $\mathbf{V}$ along $x$, $y$ and $z$ are $V_1$, $V_2$ and $V_3$.

## Step 3: what the components say

Each value is printed with its terms in post-Newtonian order, lowest first, by the printer's `rank`, so that Newton's term leads: $\Gamma^x{}_{tt} = -(\partial_xU - 4U\partial_xU + 4\partial_tV_1 + \tfrac{1}{2}\partial_t^2\partial_x\chi + \partial_x\psi)$.
The inverse is written whole: with $A$ the bracket of $g_{tt}$, $S = 1 + 2U$, $B_i = g_{ti}$ and $D = AS + B_iB_i$, $g^{tt} = -S/D$, $g^{ti} = B_i/D$ and $g^{ij} = \delta_{ij}/S - B_iB_j/(SD)$.
The Kretschmann scalar is $8\,\partial_i\partial_jU\,\partial_i\partial_jU + 4(\nabla^2U)^2$, Schwarzschild's $48m^2/r^6$ for one body.
With $U = m/r$, $\chi = mr$ and nothing moving, both charts are the published Cartesian chart of `ppn_metric` at $\beta = \gamma = 1$, which the test holds.

## Step 4: the equations of motion

The paper's (17.2) is not on the page, since it is no part of the metric, and the test file types it out to hold what the History and the drawings rest on.
In the frame of the centre of mass it is Will's (2014) relative acceleration to first post-Newtonian order, $(m/r^2)\{-\mathbf{n} + [(4 + 2\eta)m/r - (1 + 3\eta)v^2 + \tfrac{3}{2}\eta\dot r^2]\mathbf{n} + (4 - 2\eta)\dot r\,\mathbf{v}\}$.
A circular orbit closes at $\Omega^2 = (m/d^3)(1 - (3 - \eta)m/d)$ to first order in $m/d$.
An eccentric orbit integrated between two passages of its pericentre turns by $6\pi m/p$, Robertson's result, to the three percent the next order allows at $m/p = 1/400$.

## Step 5: the drawings

Every drawing takes two bodies of equal mass $m = 1$ on a circular orbit in the plane $z = 0$, a distance $d = 20\,m$ apart, at $x = \pm 10\,m$ at $t = 0$.
The angular velocity is the first order one, $\Omega^2 = (2m/d^3)(1 - 11m/2d) = 29/160000$, so each body moves at $\sqrt{29}\,c/40 = 0.135\,c$ and $\psi = (3v^2/2 - m/d)U = -73U/3200$; `nr.EIH_BINARY` holds the six potentials as functions of $t$, $x$, $y$ and $z$.
At $d = 10\,m$ the first order rate and the rate at which (17.2) itself closes a circle differ by a factor of three, so the pair is drawn no closer than this, where they differ by eight percent.
Alvi's outer limit of the near zone is $c/2\Omega = 37.1\,m$, and no drawing reaches past $30\,m$.

The spacetime diagrams draw the plane of $t$ and $z$ on the axis of the orbit, $x = y = 0$, in both charts.
Half a turn about that axis carries two equal bodies onto each other at every moment, so a ray launched along the axis stays on it.
On the axis $U = 2m/\sqrt{z^2 + 100m^2}$, $\mathbf{V} = 0$ and $\partial_t^2\chi = \partial_t\partial_z\chi = 0$, so the plane is the same in both charts, $dz/dt = \pm c\sqrt{(1 - 2U + 2U^2 - 2\psi)/(1 + 2U)}$, which is $0.70\,c$ at $z = 0$; `--verify` holds the rays to $ct \pm z_*$ with `_eih_axis`, the quadrature of the reciprocal.

The embedding diagram draws the plane $x = 0$ midway between the bodies at $t = 0$, $(1 + 2U)(dy^2 + dz^2)$ with $U = 2m/\sqrt{s^2 + 100m^2}$ at the distance $s$ from the line through the bodies, turned about that line: a surface of revolution with $\rho = s\sqrt{1 + 2U}$, flat where it crosses the line.
On that plane at $t = 0$ the two charts' times agree, since $\partial_t\chi$ vanishes there, and the generator holds the standard chart's slice to the same metric.
At another moment the plane is no longer midway between the bodies and is no surface of revolution, so one moment is drawn.

No conformal diagram is drawn: the field holds in the near zone alone, between the bodies and the distance $c/2\Omega$, which has no infinity and no horizon to place.
No figure in three dimensions is drawn: at one scale for space and time the world lines of bodies moving at $0.135\,c$ are close to vertical lines, and the turning of an orbit shows only over many revolutions.
