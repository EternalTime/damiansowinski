# Krasnikov's time machine

S. V. Krasnikov, "A singularity-free WEC-respecting time machine", Class. Quantum Grav. 15, 997 (1998), gr-qc/9702031.
Every equation below is his, by his numbering; the signature is his, $(-,+,+,+)$, and no unit is fixed.

## The two charts

The time machine chart is $L_N = Q_N \times S$, his eqs. (2), (5) and (7): $ds^2 = f^{-2}(-d\tau^2 + d\chi^2) + R^2(d\theta^2 + \sin^2\theta\,d\phi^2)$ on the square $|\tau|, |\chi| < m$ of eq. (1).
His $f$ is smooth and bounded, vanishes exactly at $(\pm h, \pm h)$ with $0 < h < m/2$, eq. (3), and repeats near the upper cut its values near the lower one, $f(\tau, \chi) = f(\tau - 2h, \chi)$, eq. (4).
The four points are removed, the square is cut along $\tau = \pm h$, $|\chi| < h$, and the upper bank of each cut is glued to the lower bank of the other.
The causality violating set, the shaded region of his Fig. 1, is the hexagon $|\tau| < h$, $|\chi| < 2h - |\tau|$, bounded by the cuts and the null lines through the four points.
$R$ is held to eq. (6), $R^{-2} \ge \max_{q_m} K$ with $K = f(f_{,\chi\chi} - f_{,\tau\tau}) + f_{,\tau}^2 - f_{,\chi}^2$.

The universe chart is $N$ on $\mathbb{R} \times S^3$, eq. (9): $ds^2 = a^2[-d\tau^2 + d\chi^2 + \rho^2(d\theta^2 + \sin^2\theta\,d\phi^2)]$ with $-\pi/2 \le \chi \le \pi/2$.
On the square $a = 1/f$ and $\rho = fR$, eq. (10a); off it $a = \hat a(\tau)$ and $\rho = \hat\rho(\chi)$, eq. (10b), with $\hat\rho = \cos\chi$ for $|\chi| > n$, $n \in (m, \pi/2)$, and $m < \pi/2$, $f = f_0 < R^{-1}\cos m$ off the square, eq. (8), for feasibility.
His word "convex" for $\hat a$ and $\hat\rho$ means bending down: his eq. (11) is $\ddot\Phi \le 0$ with $\Phi = \ln a$, $\rho''/\rho \le 0$ and $1 - \rho'^2 \ge 0$, so the page states the inequalities and never the word.
$\hat a$ vanishes at the cosmological singularities, which lie outside $L_N$; the article fixes neither their times nor $\hat a$ beyond these conditions.

Dimensions: $\tau$, $\chi$ and $\rho$ are pure numbers, $f$ is an inverse length, and $R$ and $a$ are lengths.

## What the chart checks hold

`krasnikov_time_machine.py` builds both charts with $f$, $a$ and $\rho$ free and checks them at points, with no symbolic simplification.
In the time machine chart the orthonormal Einstein tensor is $\mathrm{diag}(1/R^2, -1/R^2, -K, -K)$ and the Ricci scalar is $2K + 2/R^2$, for two unrelated smooth trial functions held in for $f$.
So the energy density is $1/R^2$, its sum with the pressure along $\chi$ vanishes, its sum with either pressure on the sphere is $1/R^2 - K$, and the weak energy condition holds exactly where $K \le R^{-2}$, which eq. (6) gives on the whole square.
The dominant energy condition would need $|K| \le R^{-2}$ as well, which eq. (6) does not bound from below, so the page does not claim it, though Krasnikov's section 3.1 does.
In the universe chart, eq. (10a) is checked to give the time machine's metric, and off the square, with $a(\tau)$ and $\rho(\chi)$, the energy density and its sums with the pressures are, times $a^2$, $3\Phi'^2 + (1 - \rho'^2)/\rho^2 - 2\rho''/\rho$, $2\Phi'^2 - 2\Phi'' - 2\rho''/\rho$ and $2\Phi'^2 - 2\Phi'' + (1 - \rho'^2)/\rho^2 - \rho''/\rho$, each non-negative under eq. (11), which is his section 3.1 in the orthonormal frame.
Where $\rho = \cos\chi$ the three pressures agree and $a^2$ times the energy density is $3(\Phi'^2 + 1)$, Friedmann's closed universe in conformal time.

## The drawings

Both planes of $\tau$ and $\chi$ are conformally flat, $f^{-2}$ and $a^2$ times the flat plane, so their light rays keep $\tau \pm \chi$ for every $f$ and every $a$.

Spacetime diagrams.
Each plane is drawn with a declared function and `any_factor`, which checks that its null directions are those of the flat plane.
The time machine's plane is the square at $m = 1/2$, $h = 1/5$, with $f = (\sin^2(\pi(\tau + h)/2h) + (\chi^2 - h^2)^2)/12R$, which vanishes in the closed square only at the four points, since $\sin^2$ vanishes again only at $\tau = \pm 3h$, outside it; obeys eq. (4) exactly, since $\sin^2$ has the period $2h$; and keeps $K$ below $0.92/R^2$, so eq. (6) holds.
It marks the two cuts and the four removed points, and its caption names the hexagon of closed timelike curves.
The universe's plane runs from pole to pole, $|\chi| \le \pi/2$, over $|\tau| \le 1$ at $n = 1$, with the declared $a = e^{\sin(\tau + 2\chi)/2}$, and marks the square, the cuts, the points and the lines $|\chi| = n$; it reads no Kretschmann scalar, since the rays need only $a$ and the drawing marks no singular edge.

Embedding diagram.
The moment $\tau = 0$ of the universe on the great sphere $\theta = \pi/2$, drawn in flat space, with $a = 1$ there as the unit of length, $R = 4/5$, $m = 1/2$ and $n = 1$.
With $f = 1$ on the line $\tau = 0$, the constant value eq. (8) gives near the square's edge, $\rho = R$ on $|\chi| \le m$, $\cos\chi$ beyond $|\chi| = n$, and between them the polynomial of degree five that matches $\rho$, $\rho'$ and $\rho''$ at both ends; it bends down, $\rho'' \le 0$, and keeps $|\rho'| \le \sin 1$, so eq. (11) holds.
Eq. (11)'s $1 - \rho'^2 \ge 0$ is the condition for the moment to lie in flat space: a cylinder of radius $R$ round the equator, where the time machine stands, between two caps of the unit sphere.
The moment is marked on both spacetime diagrams, at $\tau = 0$.

Conformal diagram, NOT_DRAWN.
$\tau$ and $\chi$ are already conformal coordinates of finite range on both planes: the time machine's square, whose four removed points Krasnikov's factor sends to infinity, and the universe's strip between its poles and its singularities.
So each spacetime diagram is the conformal diagram of its plane, and a second drawing would repeat it.

Curvature figure, not drawn.
Every scalar of both charts is a function of $f$, $a$ and $\rho$ and their derivatives, which Krasnikov fixes only by eqs. (3), (4), (6), (8), (10) and (11).
On $L_N$ the Ricci scalar is $2K + 2/R^2$, at most $4/R^2$ by eq. (6), and nothing else of it is fixed, so a figure would show a chosen $f$ and none of his spacetime.
