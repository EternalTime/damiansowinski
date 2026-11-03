# The Giddings-Strominger wormhole

Giddings and Strominger's axion wormhole is a Riemannian space of four dimensions with signature $(+,+,+,+)$ and no time in it, as the Atiyah-Hitchin space is.
Equations are named by their labels in the sources' TeX, which the arXiv copies carry, and Giddings and Strominger's by the numbers of their preprint HUTP-87/A067, read in full in the KEK scan; the published paper could not be read in full, and no equation below rests on it alone.
This note records the sources of the three charts, the notation, the checks, and what the drawings draw and why two of them are not drawn.

## Step 1: no time

No coordinate is declared a time in `DIMENSIONS`, so the checker scales nothing by $c$, and each chart's line element is the chart line element as it stands.
The Riemann tensor, the contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$ and the Weyl tensor are the same expressions in any signature.
A dot in a geodesic equation is a derivative with respect to arc length.
There is no light cone and no null infinity, so `null_rays.py` records all three charts in `NOT_DRAWN` and `conformal.py` names the space in its `NOT_DRAWN`.

## Step 2: the solution and its notation

Giddings and Strominger take the ansatz (2.7), $ds^2 = dt^2 + a^2(t)\,d\Omega_3^2$, with the field strength (2.8), $H = h(t)\,\epsilon$, where $\epsilon$ is the volume form of the 3-sphere normalised by (2.9), and $h = n/(f^2a^3)$ by (2.10), so that every 3-sphere carries the flux (2.11), $\int_{S^3}H = 2\pi^2n/f^2$.
The Einstein equation is (2.12), $(3M_p^2/16\pi)(\dot a^2/a^2 - 1/a^2) = -3n^2/(a^6f^2)$, which is $\dot a^2 = 1 - r^4/a^4$ with their $r$ the throat radius, $r^2 = 4\sqrt{\pi}\,n/(fM_p)$, stated after (2.13).
Their letters are the other way round from the page's: the page writes $a$ for the throat radius and $r$ for the radius of the 3-spheres, as Hebecker, Mikhail and Soler do in `eq:GSwormhole`, $ds^2 = (1 + C/r^4)^{-1}dr^2 + r^2d\Omega_3^2$ with $C = -n^2/(24\pi^4f^2)$, so that $a^4 = -C$.
Giddings and Strominger's (2.19), Hebecker, Mikhail and Soler's $C$ and Witten's `pinto`, $a_0 = (Gf^2m^2/3\pi^3)^{1/4}$, normalise $f$ and the Planck mass differently, so the page states only that $a^2$ is in proportion to $n$ and prints no constant.

## Step 3: the areal chart

The areal chart is Hebecker, Mikhail and Soler's `eq:GSwormhole` with $C = -a^4$, $ds^2 = dr^2/(1 - a^4/r^4) + r^2d\Omega_3^2$, one half of the wormhole, $r \ge a$, since they count the wormhole as two solutions of that form, each restricted to $r > r_0$, and the domains name the throat $r = a$ as the place the two halves meet.
Every printed value is a rational function of $r$ whose sums are powers of $r^4 - a^4$, which `quartic_pretty` names whole.
The check holds, at three points and forty digits, the Ricci scalar to $-6a^4/r^6$, the mixed Einstein tensor to $(-3, 3, 3, 3)\,a^4/r^6$, the stress of the axion, consistent with Giddings and Strominger's (2.15), $M_p^2R = -16\pi f^2H^2$, and the mixed Ricci tensor to $(-6, 0, 0, 0)\,a^4/r^6$, whose eigenvalues are non-positive with one negative: their Theorem on p. 5 asks for a negative eigenvalue somewhere, and their remark on p. 6 says that for an axion "the eigenvalues are non-positive".
The printed Kretschmann scalar is $60a^8/r^{12}$, and the Weyl tensor vanishes.

## Step 4: the proper distance chart

The proper distance chart is Giddings and Strominger's (2.7) with their $t$ written $\rho$ and their $a(t)$ written $r$, $ds^2 = d\rho^2 + r^2d\Omega_3^2$, which is Witten's `yogo`, $ds^2 = dr^2 + a^2(r)\,d\Omega^2$ with $r$ over the whole line.
It leaves $r$ a free function of $\rho$, as the radius chart of the Atiyah-Hitchin space leaves its functions free, and no component assumes (2.12), $(\partial_\rho r)^2 = 1 - a^4/r^4$, which Witten writes as `winto`, $dr = \pm\,da/\sqrt{1 - a_0^4/a^4}$, and its convention states.
Giddings and Strominger's (2.13) gives $\rho$ as a function of $r$ through elliptic integrals; the check takes $\rho$ instead by mpmath's quadrature of $1/\sqrt{1 - a^4/r^4}$ from the throat.
The check holds the Ricci scalar and the mixed Einstein and Ricci tensors to Step 3's values at three points, with $r$, $\partial_\rho r = \sqrt{1 - a^4/r^4}$ and $\partial_\rho^2 r = 2a^4/r^5$ carried along by (2.12), and holds $g_{\rho\rho}$ and $g_{\chi\chi}$ to the areal chart pulled back there.

## Step 5: the conformal chart

Loges, Shiu and Sudhir write the metric in conformal time, with the zeroth-order Einstein equation `eq:hubbleEOM`, $1 + \mathcal{H}^2 = \tilde n^2/a^4$; after the rotation $\eta \to -ir$ it reads $1 - \mathcal{H}_\mathrm{E}^2 = \tilde n^2/a^4$, with the solution $a(r) = \sqrt{\tilde n\cosh(2r)}$: "This is the GS wormhole; the wormhole throat is at $r=0$".
With their Euclidean $r$ written $\eta$ and $\tilde n = a^2$, the throat radius squared, the line element is $a^2\cosh(2\eta)(d\eta^2 + d\Omega_3^2)$, which covers both halves.
Along $r^2 = a^2\cosh 2\eta$, $r^4 - a^4 = a^4\sinh^2 2\eta$ and $dr = a^2\sinh 2\eta\,d\eta/r$, so $dr^2/(1 - a^4/r^4) = r^2d\eta^2$: the space is conformal to the cylinder $\mathbb{R} \times S^3$, so conformally flat, and its Weyl tensor vanishes in every chart.
The check holds $g_{\eta\eta}$ and $g_{\chi\chi}$ to the areal chart pulled back along that map at three points, and the Ricci scalar and the mixed Einstein and Ricci tensors to Step 3's values there.

## Step 6: the drawings

`giddings_strominger_drawings.py` draws them, and records in `NOT_DRAWN` for `null_rays.py` and `conformal.py`, as Step 1 says, that a positive definite metric has no nonzero null tangent and no causal boundary, so the space has no spacetime diagram and no conformal diagram.
The embedding diagram is the surface of $r$ and $\phi$ at $\chi = \theta = \pi/2$ of the areal chart at $a = 1$, where $g_{rr} = r^4/(r^4 - a^4)$ and $g_{\phi\phi} = r^2$, so that $\rho = r$ and $dz/dr = a^2/\sqrt{r^4 - a^4}$.
Expanding the integrand from $r$ to infinity term by term gives the height above the throat, $z = z_\infty - (a^2/r)\,{}_2F_1(\tfrac12, \tfrac14; \tfrac54; a^4/r^4)$, and Gauss's sum ${}_2F_1(\tfrac12, \tfrac14; \tfrac54; 1) = \Gamma(\tfrac54)\Gamma(\tfrac12)/\Gamma(\tfrac34)$ makes it vanish at the throat with $z_\infty = \sqrt{\pi}\,\Gamma(1/4)\,a/(4\Gamma(3/4)) = 1.3110287771\,a$.
The integrand falls as $a^2/r^2$, so the height is finite, where Ellis and Bronnikov's catenoid climbs as $\ell\,\mathrm{arccosh}(r/\ell)$ without bound.
Both halves are drawn through the throat out to $r = 4a$, the other half the same surface turned over; the generator checks $z_\infty$ against its closed form and against a quadrature in thirty digits, with $r = 1 + u^2$ taking the inverse square root at the throat out of the integrand, and each half against the closed form of $z$ and against $\rho = r$, and the two halves to meet at the throat.
The curvature figure is the stored $K$ and $R$ of the conformal chart, the one chart that covers both halves, against $\eta$ from $-1.5$ to $1.5$, where the areal radius is $a\sqrt{\cosh 3} = 3.17a$: $K = 60/(a^4\cosh^6 2\eta)$, which is $60a^8/r^{12}$, logarithmic from $10^{-5}$ to $10^2$, and $R = -6/(a^2\cosh^3 2\eta)$, which is $-6a^4/r^6$, signed logarithmic from $-10$ to $0$; both are even in $\eta$, greatest in size at the throat, and lie wholly inside their panels.
`_tools/test_giddings_strominger.py` holds the published files to the curvature scalar, both values of the Einstein tensor, the sign of the Ricci tensor's eigenvalues, the funnel's height in closed form and the drawings' registrations.
