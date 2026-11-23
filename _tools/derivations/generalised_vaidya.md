# The generalised Vaidya spacetime

Vaidya's metric with a mass function of a null time and the areal radius.
The charts are written in generalised_vaidya.py and the drawings in generalised_vaidya_drawings.py beside this file.

## The charts

### Ingoing and outgoing

Anzhong Wang and Yumei Wu, gr-qc/9803038, eq. (1), write the general spherically symmetric line element as
$ds^2 = -e^{2\psi}(1 - 2m/r)\,dv^2 + 2\epsilon e^{\psi}\,dv\,dr + r^2d\Omega^2$ with $\epsilon = \pm 1$, and take $\psi = 0$.
With $\epsilon = +1$ the null coordinate is the advanced time and the rays of constant $v$ fall in; with $\epsilon = -1$ it is the retarded time, written $u$ here.
The ingoing chart is Viqar Husain's eq. (1) with $\psi = 0$, gr-qc/9511011, and the outgoing chart is Edward N. Glass and Jean P. Krisch's eq. (1), gr-qc/9803040, $ds^2 = A\,du^2 + 2\,du\,dr - r^2d\Omega^2$ with $A = 1 - 2m/r$, turned from their signature $(+,-,-,-)$ to $(-,+,+,+)$.
Both are written with $c$, $m$ a length, and the null time a time; the components are taken in the chart $x^0 = cv$ or $x^0 = cu$.

Wang and Wu's eq. (2) is the Einstein tensor, $G^v{}_v = G^r{}_r = -2\partial_r m/r^2$, $G^r{}_v = 2\partial_v m/r^2$ and $G^\theta{}_\theta = G^\phi{}_\phi = -\partial_r^2 m/r$, whatever the sign of $\epsilon$.
The chart check holds both charts to it at three points of the plane for a mass function of both coordinates, $m = 1/2 + v/5 + r^2/9 + r\sin v/7 + vr^3/31$, which is no special case.
It holds them as well to the published charts of Vaidya's metric of the same names at $G = c = 1$ for a mass function of the null time alone, Wang and Wu's remark that $\rho = P = 0$ is Vaidya's solution, and to de Sitter space for $m = \Lambda r^3/6$, their case ii), where $G_{ab} = -\Lambda g_{ab}$.

### Husain's null fluid

Husain imposes $P = k\rho$ on the eigenvalues of the stress tensor, his eqs. (5), (7) and (11), and finds for $k \neq 1/2$ the mass function of his eq. (12), $m = f(v) - g(v)/((2k - 1)r^{2k - 1})$, and the line element of his eq. (15).
His energy density is $\rho = g/(4\pi r^{2k + 2})$ in units $G = 1$, his eq. (13), so $g \ge 0$ for positive energy, which his dominant energy condition asks for.
The chart writes his $g$ as $q^{2k}$ with $q$ a length, so that the term in the lapse is the pure number $h = (q/r)^{2k}$, as Kiselev's chart writes its term $(r_q/r)^{3w + 1}$; a length raised to the power $2k$ has no dimension the checker can balance, and $g = q^{2k} \ge 0$ is the branch the energy conditions select.
The line element is then $-\left(1 - 2f/r + 2h/(2k - 1)\right)c^2dv^2 + 2c\,dv\,dr + r^2d\Omega^2$.
The name $h$ is held as a function of $v$ and $r$ with its two first derivatives, $\partial_r h = -2kh/r$ and $\partial_v h = 2kh\,q'/q$, so every component is printed in $h$; along $v$, a time, the rate is $c$ times the derivative along $x^0$.

The chart check holds Husain's stress at $k = 3/4$ and $k = 1/3$: $G^v{}_v = G^r{}_r = -8\pi\rho$, $G^\theta{}_\theta = G^\phi{}_\phi = 8\pi k\rho$ and $G^r{}_v = 2\partial_v m/r^2$, with test functions $f = 1/2 + v/7 + v^2/40$ and $q = 1/5 + \sin^2 v/9$.
With $q = 0$ the chart is Vaidya's ingoing chart with $m = f$, and at $k = 1$ it is Bonnor and Vaidya's ingoing chart with their charge $\sqrt{2}\,q$, Husain's statement that $k = 1$ is the charged Vaidya metric.
The Kretschmann scalar is printed as Schwarzschild's $48f^2/r^6$, the cross term and the fluid's own; at $k = 1$ with $2q^2 = Q^2$ it is Reissner and Nordström's $(48M^2r^2 - 96MQ^2r + 56Q^4)/r^8$.

## The drawings

Every drawing is in units of a length $L$, and the null time enters the tanh of each source as $2cv/L$, so that the mass settles within the window drawn.

### The spacetime diagrams

The ingoing chart draws Husain's collapse of radiation, his eq. (21) at $k = 1/3$ with $C = 0$ and $A = L/2$, $B = L^{2/3}/8$: $m = (1 + \tanh(2cv/L))(L/4 + L^{2/3}r^{1/3}/16)$.
It is flat as $v \to -\infty$, and at late times its apparent horizon solves Husain's cubic $(r - 2A)^3 = (2B)^3r$, $r = 1.27079L$.
The event horizon is seeded on that radius at $cv = 2.6L$, where $\tanh(5.2)$ is $1$ to five parts in a hundred thousand, and traced back; it leaves the centre near $cv = -2.6L$.
The cones are oriented by $cv - r$, a time wherever $g^{rr} = 1 - 2m/r < 2$, which holds since $m \ge 0$.

Husain's chart draws his eq. (20), $k = 3/4$, $A = L$ and $B = 1/2$: $f = L(1 + \tanh(2cv/L))/2$ and $q = L((1 - \tanh(2cv/L)/2)/4)^{2/3}$, so that $2f/r$ and $4h = (1 - B\tanh)/r^{3/2}$ are his two terms.
At early times the lapse is $1 + (3/2)(L/r)^{3/2}$ and the centre is a naked singularity; at late times it is $1 - 2L/r + (L/r)^{3/2}/2$, a black hole whose horizons are $r = 1.60538L$ and $0.06690L$.
Near the early singularity $g^{rr}$ exceeds $2$ and $cv - r$ is no time, so that view orients its cones by the ingoing rays, future directed in every ingoing chart.

The outgoing chart draws Glass and Krisch's diffusing strings, their eq. (18a) of the letter and (38) of the paper, with $\rho_0 = 0$: $m = m_0 + 2\pi k_1(r^2 - 2Dcu)$ at $m_0 = L/2$, $k_1 = 1/80\pi L$ and $D = L/2$.
The strings' density is $k_1/r$, their eq. (17a), and the mass falls by $4\pi Dk_1$ per unit of $cu$, so the star loses its mass to the light while the strings diffuse inward.
Inside $r = 2m$ the chart covers a white hole, as Vaidya's outgoing chart does.

### No conformal diagram

The mass function is free in both coordinates, and a conformal compactification needs the global structure of one member: its event horizons and its null infinity are fixed by the whole history of $m$, as the Mallett radiating star's are by its free mass history, which is drawn without one too.
No source gives a member of the family a compactification in closed form.

### The embedding diagram

The moments are the slices of constant $cv - r = TL$ of the collapse of radiation, spacelike since $g^{rr} < 2$; with $v = T + r$ the metric pulls back to $(1 + 2m/r)\,dr^2 + r^2d\phi^2$, so $dz/dr = \sqrt{2m/r}$.
At every finite $v$ the centre holds the mass $m(v, 0) = L(1 + \tanh(2cv/L))/4$, and the surface closes there in a spike.
The marginally trapped sphere on each moment is the root of $r = 2m(T + r, r)$, found by bisection and checked.
The movie runs from $T = -1$ to $1.5$ in steps of $1/16$, forty one frames; before $T = -1$ the central mass is below a hundredth of $L$ and its spike narrower than the profile is sampled, where the isometry check fails along the first chord.
The moments lie on the ingoing chart's spacetime diagram as the level lines of its vertical axis; the other two diagrams draw other spacetimes, and the exceptions file names them.

### The curvature figure

The published Kretschmann and Ricci scalars of the ingoing chart over the plane of $r$ and $v$, from $r = 0$ to $3L$ and $cv = -1.5L$ to $1.5L$, with the collapse of radiation declared, and the marginally trapped spheres where the published $g^{rr}$ vanishes.
The Ricci scalar there is $(4\partial_r m + 2r\,\partial_r^2 m)/r^2 = (1 + \tanh(2cv/L))\,L^{2/3}r^{-8/3}/18$, the trace of the fluid, which the null dust alone does not have; the figure checks it at every cell.
