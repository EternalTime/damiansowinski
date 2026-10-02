# The Neugebauer-Meinel disc

Neugebauer and Meinel's rigidly rotating disc of dust, Phys. Rev. Lett. 75, 3046 (1995), arXiv:gr-qc/0302060, in five charts.
`_tools/derivations/print_charts.py --metric neugebauer_meinel` writes them in about three minutes, and `verify_metrics.py --system neugebauer_meinel/<chart>` checks each, the spheroidal chart in two minutes and the others in seconds.
`_tools/derivations/nm_disc.py` evaluates the solution itself, and its opening lines are the full statement of how.

## Step 1. The charts and their sources

`weyl` is the metric as the letter of 1995 writes it, its (1), in the coordinates of Weyl, Lewis, and Papapetrou: $ds^2 = e^{-2U}(e^{2k}(d\rho^2 + dz^2) + \rho^2d\phi^2) - e^{2U}(c\,dt + a\,d\phi)^2$.
Their $\zeta$ is written $z$, since the checker's reader takes no derivative along a coordinate named $\zeta$.
`corotating` is the same form in the frame that turns with the disc, $\varphi = \phi - \Omega t$, with their $U'$, $a'$ and $k'$, the footnote of Meinel's lecture of 1997 (arXiv:gr-qc/9703077) and section 2 of Neugebauer, Kleinwächter and Meinel, Helv. Phys. Acta 69, 472 (1996), arXiv:gr-qc/0301107.
`bardeen_wagoner` is Bardeen and Wagoner's form, Astrophys. J. 158, L65 (1969), their (1), and 167, 359 (1971), their (II.1) with $B = 1$, a lapse $e^\nu$ and a rate of dragging $\omega$; their $\mu$ is written $\alpha$, as Meinel writes it in Ann. Phys. 11, 509 (2002), since $\mu$ is the parameter of the family.
`spheroidal` is that form in the oblate spheroidal coordinates both pairs of authors compute in, $\rho = \rho_0\sqrt{(1 + \xi^2)(1 - \eta^2)}$ and $z = \rho_0\xi\eta$, Bardeen and Wagoner's letter and the letter of 1995, its (11).
`black_hole_limit` is the field outside the disc at $\mu \to \mu_0$, the extreme Kerr metric as Bardeen and Wagoner's (VIII.2) to (VIII.5) of 1971 write it, with $r + m$ the radius of Boyer and Lindquist; Meinel's (33) of 2002 is its Ernst potential.

## Step 2. Why the functions are left free

The solution's $U$, $a$ and $k$ are quotients of theta functions of a curve of genus two whose branch points move with $\rho$ and $z$, and no component written out in them could be printed or read.
So the four charts of the disc leave their three functions free, as Gowdy's, Black Saturn's and the double Kerr solution's charts do, no component assumes a field equation, and the parameters' descriptions state the vacuum equations that hold off the disc.
With $x^0 = ct$ the function $a$ is a length and $\omega$ an inverse length, the dragging per unit of $ct$.
`OVERRULED` in `metric_tags.py` sets nothing: the disc is dust, and the spacetime is no vacuum.
Inside the ergoregion, which the disc has for $\mu > 1.68849$, $e^{2U}$ is negative and $U$, $a$ and $k$ are not real, which is why the drawings of the plane read Bardeen and Wagoner's form, whose three functions are real there.

## Step 3. What the script checks

`neugebauer_meinel_check` holds the equations each free chart states to making every Ricci component vanish and to an integrable quadrature: Ernst's equation written in $U$ and $a$ with the quadrature for $k$, and Bardeen and Wagoner's two equations with the quadrature for $\alpha$, which is
$$\partial_\rho\alpha = -\partial_\rho\nu + \rho\left((\partial_\rho\nu)^2 - (\partial_z\nu)^2\right) - \tfrac{1}{4}\rho^3e^{-4\nu}\left((\partial_\rho\omega)^2 - (\partial_z\omega)^2\right), \qquad \partial_z\alpha = -\partial_z\nu + 2\rho\,\partial_\rho\nu\,\partial_z\nu - \tfrac{1}{2}\rho^3e^{-4\nu}\partial_\rho\omega\,\partial_z\omega.$$
The corotating chart is held to being Weyl's with $\phi = \varphi + \Omega t$ and $e^{2U'} = e^{2U}((1 + \Omega a/c)^2 - \Omega^2\rho^2e^{-4U}/c^2)$, $(1 - \Omega a'/c)e^{2U'} = (1 + \Omega a/c)e^{2U}$, $k' - U' = k - U$.
Bardeen and Wagoner's chart is held to being Weyl's with $e^{2\nu} = \rho^2e^{2U}/(\rho^2 - a^2e^{4U})$, $\omega = a\,e^{4U}/(\rho^2 - a^2e^{4U})$ and $\alpha = k - U$, and the spheroidal chart to being Bardeen and Wagoner's pulled back.
The limit chart is held to a vanishing Ricci tensor and to being the published Boyer-Lindquist chart of `kerr` at $a = GM/c^2 = m$ with $r \to r + m$.
The theta functions are held to Ernst's equation by fourth order differences at six points, and their $a$ to $\partial_\rho a = \rho\,e^{-4U}\partial_zb$.

## Step 4. The solution as numbers

`nm_disc.Point` evaluates the Ernst potential from the theta formula of Neugebauer, Kleinwächter and Meinel's (2.44), which the review of 2003 repeats as its (113): Gauss-Legendre quadrature along polylines between the six branch points, the square root followed continuously along each.
Their figure of the Riemann surface shows the cut of $z \mp i\rho$ to the right of the cuts of $X_1$ and $X_2$; nearer the disc that pair passes over and under the cut of $X_2$, or between its ends, and `Point.ways` holds the curves of each arrangement.
The metric is taken from the appendix of the review of 2003, J. Math. Phys. 44, 3407.
As printed, its second formula has $\vartheta^*(\mathbf{0})$ in the denominator where $\vartheta^*(\mathbf{c})$ belongs: with that put right the quotient and its companion, with $\vartheta$ and $\vartheta^*$ exchanged in two places, sum to 2 to ten digits at every point tried, and agree with a quadrature of $\partial_\rho a = \rho\,e^{-4U}\partial_zb$.
The three functions are written in combinations that stay finite on the ergosurface: $F = e^{2U}$, $A = a\,e^{2U} = -g_{t\phi}$, $g_{\phi\phi}$, and $g_{\rho\rho} = e^{2k - 2U}$.
The theta functions are never told the conditions on the disc, and satisfy them: $e^{2U'} = e^{2V_0}$ at every radius of the disc to nine digits, the metric of the disc worked out from the two conditions alone agrees with theirs to seven, $V_0$ is the closed form of 1994 with Weierstrass's function, and $e^{2U} = 1 - \mu/2$ at the rim.
`_tools/test_neugebauer_meinel.py` holds all of it, with the Maclaurin disc at $\mu = 0.01$ and the extreme Kerr potential near $\mu_0$.

## Step 5. The disc that is drawn

Every drawing of the disc is $\mu = 3$, the disc whose Ernst potential Neugebauer, Kleinwächter and Meinel plot, in units of its coordinate radius $\rho_0$.
There $e^{2V_0} = 0.0303$, a redshift of $4.74$ from the centre, $\Omega\rho_0 = 0.2133\,c$, $GM/c^2 = 1.716\,\rho_0$, the ergosurface crosses the plane at $0.150$ and $1.831\,\rho_0$, and a point at rest in the turning frame moves at the speed of light at $1.378\,\rho_0$.
`nm_disc.plane` tabulates $F$, $A$, $g_{\phi\phi}/\rho^2$ and $g_{\rho\rho}$ on the plane at Chebyshev points, in $\eta = \sqrt{1 - \rho^2/\rho_0^2}$ on the disc and in $\xi/(1 + \xi)$ outside it, since the functions are smooth in the spheroidal coordinates and have a term in $(\rho - \rho_0)^{3/2}$ at the rim in Weyl's.
`null_rays.DECLARED_FUNCTIONS` carries the tables into the generators as sympy functions with their derivatives.

## Step 6. The diagrams

The rotations fix the axis, and the reflection $z \to -z$ fixes the plane of the disc, so both are totally geodesic.
The spacetime diagrams draw the axis in Weyl's form, where $e^{2U}$ is positive; the plane in Bardeen and Wagoner's form with $\phi$ divided out, across the ergoregion; the plane at rest in the turning frame out to the radius where it moves at the speed of light; the axis and the plane beyond the rim in the spheroidal chart; and the axis and the equator of the limit.
`--verify` checks the rays against tortoise coordinates taken by quadrature, and the limit's axis against $r_* = r + 2m\ln(r/m) - 2m^2/r$.
The embedding diagram is the plane of the disc at one moment, circles of radius $\rho\,e^{-\nu}$.
The circles reach $4.41\,\rho_0$ at $0.92\,\rho_0$, narrow to $4.07\,\rho_0$ at $1.20\,\rho_0$ and widen again, the beginning of the throat of the limit.
Where they change faster than the distance across them the slice stands in Minkowski space: inside $\rho = 0.689\,\rho_0$ and between $0.979$ and $1.006\,\rho_0$; the four pieces join on three level circles.
The conformal diagram draws the whole axis through the disc, a diamond with null infinity on every side, and the axis of the limit, a diamond with the horizon on its left.
The limit's views mark no moment of the embedding, since the limit is another spacetime than the disc drawn.
