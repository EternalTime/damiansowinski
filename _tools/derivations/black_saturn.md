# Black Saturn

Elvang and Figueras's black Saturn, JHEP 05 (2007) 050, arXiv:hep-th/0701035: a black ring in balance around a spherical black hole in five dimensions.
Three charts are published, Weyl's canonical coordinates, the polar coordinates of their section 3.3, and the ring coordinates of their appendix A.2 for the ring alone.

## Step 1: the units and the parameters

Their Weyl coordinates $\rho$ and $z$ are areas, and their section 3.1 takes the scale $L^2 = a_2 - a_1$ out of them.
The charts use $\rho/L^2$ and $\bar z = (z - a_1)/L^2$, pure numbers, written $\rho$ and $z$, so the rod ends $a_1, a_5, a_4, a_3, a_2$ are $0, \kappa_3, \kappa_2, \kappa_1, 1$ and the one length $L$ carries every dimension.
$\beta$ is their $\bar c_2$, renamed because the checker's reader takes no accent on a parameter.
In units of $L$, $c_1 = \sqrt{2\kappa_1\kappa_2/\kappa_3}$ is their condition for a regular axis, $c_2 = \beta c_1(1 - \kappa_2)$ their definition of $\bar c_2$, $q = \beta c_1/(1 + \kappa_2\beta)$ their $q$, and $k = 1/(1 + \kappa_2\beta)$ their $k$ without its absolute value, since only $k^2$ enters.
The solitons $\mu_i$, the polynomials $M_0$ to $M_4$ and $F$, and $H_x$, $H_y$, $G_y$, $P$ and $\omega_\psi$ are those of their section 2.3, each a name the charts define.
$1/\sqrt{G_x}$ in their $\omega_\psi$ is written $\sqrt{G_y/\rho^2}$, since $G_xG_y = \rho^2$.

## Step 2: the form of the metric

Their own form of the metric, in section 2.3, is written around $H_y/H_x = -g_{tt}$, which changes sign on the ergosurfaces, so a chart built on $e^{2U} = H_y/H_x$ is complex inside every ergoregion, and one built on their own six functions takes the chart printer more than a quarter of an hour.
The charts take the same metric as a lapse and a shift,
$$ds^2 = -e^{2W - 2V}c^2dt^2 + L^2\left(e^{2V}\left(d\psi - \frac{\Omega}{L}c\,dt\right)^2 + \rho^2e^{-2W}d\phi^2 + e^{2\nu}\left(d\rho^2 + dz^2\right)\right),$$
which is the general form of their section 2.1 with $\det G = -\rho^2$ and is real wherever the circles of $\psi$ are spacelike.
$V$, $\Omega$, $W$ and $\nu$ are free functions in every component, as Gowdy's three are, and for black Saturn
$$e^{2V} = \frac{G_yH_x^2 - (\omega_\psi + qH_y)^2}{H_xH_y}, \qquad \Omega = \frac{\omega_\psi + qH_y}{H_x}e^{-2V}, \qquad e^{2W} = G_y, \qquad e^{2\nu} = k^2H_xP.$$
The vacuum equations are Laplace's equation for $W$, $\nabla^2V = -\tfrac12e^{4V - 2W}|\nabla\Omega|^2$, $\nabla^2\Omega = -(4\nabla V - 2\nabla W)\cdot\nabla\Omega$, and a quadrature for $\nu$, with $\nabla^2 = \partial_\rho^2 + \rho^{-1}\partial_\rho + \partial_z^2$.

## Step 3: what the script checks

`black_saturn_check` in `print_charts.py` holds Weyl's chart to three things.
The stated equations make every Ricci component vanish, symbolically.
Elvang and Figueras's functions, read from the names the chart defines and so from the text the page prints, make every Ricci component vanish at three points in forty digits, for a Saturn out of balance with a spinning hole, $\kappa = (7/10, 9/20, 1/5)$ and $\beta = 3/10$; they write that they checked the vacuum equations numerically themselves.
At $r = 1000L$ the four functions are flat space's.
The polar chart, $\rho = r^2\sin\theta\cos\theta/L^2$ and $z = r^2(\cos^2\theta - \sin^2\theta)/(2L^2)$, is held to being Weyl's pulled back.

## Step 4: the ring alone

The ring chart is the ring metric of their appendix A.2 with $\psi \to -\psi$, the sense of rotation they give the Saturn itself, so that all three charts turn toward increasing $\psi$.
Its Ricci tensor vanishes exactly.
It is Weyl's chart at $\kappa_1 = 1$ and $\beta = 0$ through the map of that appendix: with $\alpha = (\nu(1 + \lambda) - 2\lambda)/(2(1 - \lambda))$ the rod ends in units of $R^2$ are $\alpha$, $-\nu/2$, $\nu/2$ and $1/2$, so $L_E^2 = (1/2 - \alpha)R^2$.
Their appendix takes $k^2 = (1 - \lambda)/(1 - \nu)^2$ where the Saturn has $k = 1$, which rescales the angles: $\psi$ and $\phi$ of the ring chart have the period $2\pi k$, and the Saturn's scale is $L = kL_E$.
The script compares the two metrics slot by slot at three points for a ring out of balance, to twenty digits.

## Step 5: the plane of the ring

The drawings take a hole with no angular momentum of its own inside a ring in balance: $\kappa = (7/8, 9/16, 3/7)$ and $\beta = 0$, where $(\kappa_1 - \kappa_2)^2 = \kappa_1(1 - \kappa_2)(1 - \kappa_3)(\kappa_1 - \kappa_3)$ and $c_1^2 = 147/64$.
On $\rho = 0$ each soliton is $2(a_i - z) + O(\rho^2)$ below its rod end and $\rho^2/(2(z - a_i)) + O(\rho^4)$ above it, and the leading terms of the four functions are rational in $z$.
With $Q = 939 - 2527z + 2345z^2 - 784z^3$, outside the ring, $z < 3/7$,
$$e^{2V} = \frac{2Q}{49(1 - z)(9 - 16z)}, \quad \Omega = \frac{21\sqrt3\,(7 - 8z)}{Q}, \quad e^{2W} = \frac{4(3 - 7z)(7 - 8z)}{7(9 - 16z)}, \quad e^{2\nu} = \frac{7(9 - 16z)}{4(3 - 7z)(7 - 8z)},$$
and between the ring and the hole, $9/16 < z < 7/8$,
$$e^{2V} = \frac{14z(1 - z)}{7z - 3}, \quad \Omega = \frac{3\sqrt3}{16z}, \quad e^{2W} = \frac{7(7 - 8z)(16z - 9)}{64(7z - 3)}, \quad e^{2\nu} = \frac{64(7z - 3)}{7(7 - 8z)(16z - 9)}.$$
$\Omega$ is $1/\sqrt3$ at both edges of the ring and $3\sqrt3/14$ at the hole, their horizon angular velocities at these values.
$g_{tt}$ vanishes at $z = 0$ outside the ring and is positive across the whole gap, so the two ergoregions have merged, as in Elvang, Figueras, Horowitz, Hubeny and Rangamani, Class. Quantum Grav. 26, 085011 (2009).
`_tools/test_black_saturn.py` holds these forms to the functions the page prints, at $\rho = 10^{-12}$ in sixty digits.

## Step 6: the diagrams

The plane of the ring is where the circles of $\phi$ have shrunk to points, and the circles of $\psi$ are divided out of it as for the BTZ hole: the rays of no angular momentum run at $dz/d(ct) = \pm e^{W - V - \nu}/L$.
The rows read no Kretschmann scalar, whose second derivatives across the plane the plane's own functions do not hold.
`quotient_checks` counts a component that is not finite on an axis and multiplies a vanishing part of the tangent as nothing, since $\Gamma^\phi{}_{\rho\phi} = 1/\rho - \partial_\rho W$ there.
The plane across the ring's, $\rho = 0$ and $z > 1$, has no view: the circles of $\psi$ have shrunk there, where the lapse and shift form ends.
The embedding diagram is the same plane at constant $t$, the surface of $z$ and $\psi$ with $g_{zz} = L^2e^{2\nu}$ and circles of radius $Le^V$.
Near the ring $g_{zz} - (d\rho/dz)^2$ is negative, from $z = 0.242$ to within $0.005$ of the outer edge and from within $0.005$ of the inner edge to $z = 0.787$, so flat space carries the plane only outside the first circle and the band about the hole inside the second, which ends on the circle of radius $7L/10$ where the hole's horizon meets the plane.
The ring's two edges have the radii $\sqrt{15/2}\,L$ and $\sqrt{147/40}\,L$.
No conformal diagram is drawn, since no published chart covers the inside of either horizon.
The ring alone is drawn at $\nu = 1/2$, $\lambda = 4/5$, the balanced ring of least angular momentum for its mass, $j^2 = 27/32$.
