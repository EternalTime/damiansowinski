# The point charge in three dimensions

The static spacetime of a point electric charge in two dimensions of space and one of time, with signature $(-,+,+)$.
This note records the sources of its three charts, how each author's constants are carried onto the page's $k$, $Q$ and $r_c$, the chart that is left out and why, and what the drawings draw.
Gott, Simon and Alpert's paper of 1986, Deser and Mazur's of 1985 and Melvin's of 1986 could not be read in full; their abstracts, the zbMATH review of the first and the papers below that restate the solution were read, and no equation here rests on an unread text alone.

## Step 1: the constants

The line element of the static chart is Gott, Simon and Alpert's as their abstract prints it, $ds^2 = -(kQ^2/2\pi)\ln(r_c/r)dt^2 + (2\pi/kQ^2)[\ln(r_c/r)]^{-1}dr^2 + r^2d\phi^2$, with the field $E = Q/r$ and $r_c$ "a constant"; I. Gottlieb's review in zbMATH (Zbl 0604.53037) calls $r_c$ "an integration constant, independent of Q".
Cataldo and Salgado print the same metric in signature $(+,-,-)$ (their equation `ec.2`) and take from Gott, Simon and Alpert the stress tensor $T_{(a)(b)} = g_{(a)(b)}F^2/8\pi - F_{(a)(c)}F_{(b)}{}^{(c)}/2\pi$ (their equation `tensor electromagnetico`), which in signature $(-,+,+)$ is $T_{\mu\nu} = (F_{\mu\alpha}F_\nu{}^\alpha - g_{\mu\nu}F^2/4)/2\pi$.
Their point mass, $ds^2 = dt^2 - dr^2/(1 - kM/2\pi)^2 - r^2d\phi^2$ (their equation `ec.1`), cited to the same paper, has the angle deficit $kM$, so $k$ is the constant of $G_{\mu\nu} = kT_{\mu\nu}$.
The page writes $k = 8\pi G/c^4$ with $G$ Newton's constant in three dimensions, which gives a particle of mass $m$ and rest energy $mc^2$ the deficit $8\pi Gm/c^2$ of the point particle's page.
With $F = (Q/r)\,dr \wedge c\,dt$ the metric solves $G_{\mu\nu} = kT_{\mu\nu}$ and Maxwell's equations exactly, which each chart's check confirms at points, and the field in an orthonormal frame is $E = Q/r$.
$kQ^2$ is a number: $k$ carries $\mathrm{T}^2/(\mathrm{M}\,\mathrm{L}^2)$ and $Q$ carries $\mathrm{M}^{1/2}\mathrm{L}/\mathrm{T}$, as `register_dimensions` declares.
The radial component is printed as $2\pi/(kQ^2\ln(r_c/r))$, the same as the abstract's $[\ln(r_c/r)]^{-1}$ factor, since the reader takes no square bracket as a group.

## Step 2: the outgoing chart

Podolský and Papajčík's charged black holes of the Robinson-Trautman class are $ds^2 = r^2(d\varphi + \alpha\,du)^2 - 2\,du\,dr + (m - \kappa_0Q^2\ln|Q/r| + \Lambda r^2)du^2$ (their `RTmetric-final-Qconst`), with $F = (Q/r)\,dr \wedge du$ (their `RTformF-tr-Qconst`), Einstein's equations $G_{ab} + \Lambda g_{ab} = 8\pi T_{ab}$ and $T_{ab} = (\kappa_0/8\pi)(F_{ac}F_b{}^c - g_{ab}F^2/4)$.
At $\alpha = 0$ and $\Lambda = 0$ it is $-f\,du^2 - 2\,du\,dr + r^2d\varphi^2$ with $f = -m + \kappa_0Q^2\ln|Q/r|$, which is the page's $f = (kQ^2/2\pi)\ln(r_c/r)$ under $\kappa_0 = k/2\pi$ and $r_c = |Q|e^{-m/\kappa_0Q^2}$.
Their transformation $du = dt + (m - \kappa_0Q^2\ln|Q/r| + \Lambda r^2)^{-1}dr$ (their `du->dt`) is $c\,du = c\,dt - dr/f$, which the chart's check holds by pulling the static metric back at four points.
The constant of $u$ is fixed so that $u = t - r_*$ with the tortoise coordinate of Step 4.

## Step 3: the ingoing chart

Mao and Zhao write the charged rotating solution in the advanced time, $ds^2 = (r^2f^2/R^2)dv^2 - 2(r/R)dv\,dr - R^2(d\phi - N^\varphi dv)^2$ in signature $(+,-,-)$, through $dv = dT + R\,dr/(rf^2)$ (their section on the charged rotating solution in the gauge of Newman and Unti).
At $\omega = 0$ and no cosmological constant, $R = r$ and $N^\varphi = 0$, and in signature $(-,+,+)$ it is $-f\,dv^2 + 2\,dv\,dr + r^2d\phi^2$, the time reverse of Step 2, reached by $c\,dv = c\,dt + dr/f$; its $f^2 = -M - Q^2\log r^2$ is the page's $f$ with their units of charge.
The constant of $v$ is fixed so that $v = t + r_*$.

## Step 4: the chart left out

Gott, Simon and Alpert's abstract says they produced a Kruskal extension, with two static regions $r < r_c$ and two dynamical ones $r > r_c$; its coordinates are in the unread body of the paper, so no Kruskal chart is printed.
The conformal diagram builds its own map from the published static and null charts, below, and it is not a chart of the page.

## Step 5: the tortoise coordinate and the horizon

$f = a\ln(r_c/r)$ with $a = kQ^2/2\pi$ vanishes only at $r = r_c$, where $f'(r_c) = -a/r_c$ and the surface gravity is $\kappa = a/2r_c$.
With $x = \ln(r_c/r)$ inside the horizon and $y = \ln(r/r_c)$ outside it, $r_* = (r_c/a)E_1(x)$ and $r_* = -(r_c/a)\mathrm{Ei}(y)$ have $dr_*/dr = 1/f$; the first vanishes at $r = 0$, both tend to $+\infty$ at $r_c$, and the second runs to $-\infty$ as $r \to \infty$, since $\mathrm{Ei}(y) \sim e^y/y$.
So $r = 0$ lies at a finite tortoise coordinate and is a timelike line, and $r \to \infty$ lies at an infinite one, which null rays reach at an infinite affine parameter: $r$ is affine along the radial null rays of both null charts.
With $U = \pm e^{\kappa u}$ and $V = \pm e^{-\kappa v}$, $UV = F(r)$, $F = -e^{-E_1(x)}$ inside and $e^{\mathrm{Ei}(y)}$ outside, and both are $e^\gamma\ln(r/r_c)$ to first order, $\gamma$ Euler's constant, so $F$ is smooth through $r_c$ and the extension joins the regions there.
The Kretschmann scalar $K = 3k^2Q^4/(4\pi^2r^4)$ and the Ricci scalar $R = kQ^2/(2\pi r^2)$ are finite at $r_c$ and diverge at $r = 0$; the Cotton tensor does not vanish, so the spacetime is not conformally flat.

## Step 6: the drawings

Every drawing is made at $kQ^2 = \pi$, so $a = 1/2$, $f = \tfrac12\ln(r_c/r)$ and $\kappa = 1/4r_c$, with $r_c$ the unit of length: the slice $t = 0$ then embeds in flat space from $r = e^{-2}r_c$ to $r_c$, and the two static regions together are a closed surface whose width over height is $0.47$.
The spacetime diagrams draw the plane $\phi = 0$ of each chart, on which every radial null geodesic stays; the static chart and the outgoing chart are oriented by the outgoing rays, which reads the region beyond $r_c$ as the expanding one, and the ingoing chart by the ingoing rays, which reads it as the contracting one, the past of the static region.
`--verify` checks the rays against $t \pm r_*$, $u + 2r_*$ and $v - 2r_*$ with the closed forms of Step 5.
The conformal diagram is the whole spacetime, each point a circle, by $p = \arctan U$ and $q = \arctan V$ in the four cells of `charged_three_dimensions_drawings.py`: $r = 0$ is $UV = -1$, the lines $X = \pm\pi/2$, and $r \to \infty$ is the pair of null lines above and below, so the spacetime is a hexagon.
Each map is checked null, future directed and the inverse of the published metric at random points, against $UV = F(r)$, and against the limits $r \to 0$ and $r \to \infty$, and the published Kretschmann scalar is checked to diverge at $r = 0$ and to stay finite on the horizon.
The embedding diagram is the slice $t = 0$ of both static regions, joined at the horizon $r = r_c$, the widest circle, as the slice inside the inner horizon of the Reissner-Nordström black hole is drawn; nearer the charge than $e^{-2}r_c$ the circles grow faster than the distance out to them, and beyond the horizon $r$ is a time.
The curvature figure draws the published $K$ and $R$ of the static chart against $r$ from $0.1$ to $3\,r_c$, with the horizon marked, and checks $K = 3R^2$, $R = 1/(2r^2)$ and the fall of both outward.
