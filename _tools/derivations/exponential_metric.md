# The exponential metric of Papapetrou and Yilmaz

The static, spherically symmetric metric with $g_{tt} = -1/g_{rr} = -e^{-2m/r}$ in isotropic coordinates.
This note records the four charts, where each comes from, what each is held to, and what each diagram draws.
The sources have $G = c = 1$; the metric file keeps $c$ and writes the mass as the length $m = GM/c^2$.

## Step 1: the papers

Three papers were verified on Crossref before anything was built.
A. Papapetrou, "Eine Theorie des Gravitationsfeldes mit einer Feldfunktion", Zeitschrift für Physik 139, 518 (1954), doi:10.1007/BF01374560.
Huseyin Yilmaz, "New Approach to General Relativity", Physical Review 111, 1417 (1958), doi:10.1103/PhysRev.111.1417.
Petarpa Boonserm, Tritos Ngampitipan, Alex Simpson and Matt Visser, "Exponential metric represents a traversable wormhole", Physical Review D 98, 084048 (2018), doi:10.1103/PhysRevD.98.084048, arXiv:1805.03781.
The first two are not open, so what each did is taken from Makukov and Mychelkin, Physical Review D 98, 064050 (2018), their section II, and from the abstract of Yilmaz's paper on the publisher's page, which also gives his affiliation.
Every equation number below is of the 2018 paper of Boonserm and his coauthors unless another is named.

## Step 2: the isotropic chart

Their (1.1): $ds^2 = -e^{-2m/r}dt^2 + e^{2m/r}\left(dr^2 + r^2d\Omega^2\right)$.
The area of the sphere $r$ is $4\pi r^2e^{2m/r}$, their (2.1), with its least value $4\pi e^2m^2$ at $r = m$, the throat, their (2.2) and (2.3).
`print_charts.exponential_metric_check` holds the chart to that, to the Ricci scalar $-2m^2e^{-2m/r}/r^4$ of their (4.7), and the printed Kretschmann scalar is their (4.9).
The printed Riemann, Weyl, Ricci and Einstein tensors are their (4.1) to (4.8) with the indices placed as the collection places them.

## Step 3: the matter

Their (9.1): $R_{ab} = -\tfrac{1}{2}\nabla_a\Phi\nabla_b\Phi$ with $\Phi = 2m/r$, which is $R_{ab} = -2\,\partial_a\varphi\,\partial_b\varphi$ for $\varphi = m/r$, Einstein's equation for a massless scalar field of negative kinetic energy.
The check holds every chart to that in every slot and to the wave equation $\Box\varphi = 0$, written through the logarithmic derivative of the determinant so that no root is taken.
Fisher, Janis, Newman and Winicour's metric has $R_{ab} = +2\,\partial_a\varphi\,\partial_b\varphi$, as `fisher_jnw.md` Step 2 records, so the two differ by the sign of the field's energy.
Makukov and Mychelkin's (7) to (10) reach the exponential metric as the limit $\gamma \to \infty$ of that metric with $\gamma b = 2m$ fixed, $(1 - 2m/\gamma r)^\gamma \to e^{-2m/r}$, which is the relation the two pages state.

## Step 4: the Cartesian chart

$ds^2 = -e^{-2m/r}c^2dt^2 + e^{2m/r}\left(dx^2 + dy^2 + dz^2\right)$ with $r = \sqrt{x^2 + y^2 + z^2}$ a name the chart defines.
It is Misner's (3.1), the metric he attributes to Yilmaz, $ds^2 = -e^{2\Phi}d(ct)^2 + e^{-2\Phi}(dx^2 + dy^2 + dz^2)$, at $\Phi = -m/r$, and the form of Boonserm and his coauthors' (6.3) with the conformal factor put back.
The check pulls the isotropic chart back along $r = \sqrt{x^2 + y^2 + z^2}$, $\theta = \arccos(z/r)$, $\phi = \mathrm{atan2}(y, x)$.

## Step 5: the areal chart

Their (3.8) and (3.9): with $R = re^{m/r}$, $ds^2 = -e^{-2m/r}dt^2 + dR^2/(1 - m/r)^2 + R^2d\Omega^2$, and their (3.13) and (3.14) invert it, $r = -m/\mathrm{W}(-m/R)$.
They write $r_s$ for the areal radius, which the collection keeps for the Schwarzschild radius, so the chart writes $R$, as Makukov and Mychelkin's (25) to (27) of 2020 do.
The principal branch is the near side, $r > m$, and $\mathrm{W}_{-1}$ the far side; the areal radius has its least value $e\,m$ on the throat, so the chart covers one side at a time.
The isotropic radius is held as a function $r(R)$, `HELD` in `verify_metrics.py`, and `exponential_areal` writes every derivative of it by $\partial_R r = r^2/(R(r - m))$ and every exponential by $e^{m/r} = R/r$, so that a value is a rational function of $r$, $R$ and $m$ with no relation left among them.
The check holds that rate and $re^{m/r} = R$ to the definition, which the checker's `norm` reduces by $e^{k\mathrm{W}(x)} = (x/\mathrm{W}(x))^k$.
The line element and the two metric components it names keep the exponential, as their (3.9) writes them; every other value is printed in $r$ and $R$, so the Kretschmann scalar reads $4m^2(12r^2 - 16mr + 7m^2)/r^4R^4$.

## Step 6: the harmonic chart

Bronnikov, Fabris and Zhidenko (2011), their (53) and (54) with $s(k, v) = v$, the case $k = 0$: $ds^2 = e^{-2mv}dt^2 - (e^{2mv}/v^2)\left(dv^2/v^2 + d\Omega^2\right)$, and their (60) is the isotropic chart with $u = 1/v$.
The collection's chart of Fisher, Janis, Newman and Winicour's metric calls the harmonic coordinate $u$, and this chart does the same, so its $u$ is their $v$ and is $1/r$.
The check holds $\Box u = 0$, and the field is $\varphi = mu$.
Their Branch B is where the two statements about the ends come from: the throat "has the size $e \cdot m$", and $r \to 0$ is "a singular horizon", with all curvature invariants going to zero and no continuation.

## Step 7: the singular horizon

On the plane of $t$ and $r$, $g_{tt}g_{rr} = -1$, so along a radial null geodesic $e^{-2m/r}\,d(ct)/d\lambda = E$ is constant and $dr/d\lambda = \pm E$: $r$ is an affine parameter, and a ray reaches $r = 0$ at a finite affine distance.
The coordinate time it takes is infinite, since the tortoise coordinate $r_* = re^{2m/r} - 2m\,\mathrm{Ei}(2m/r)$ falls as $-r^2e^{2m/r}/2m$.
The Kretschmann scalar goes to zero there, but $R_{ab}k^ak^b = R_{rr}E^2 = -2m^2E^2/r^4$ along the ray does not, so the curvature in a frame carried along the ray diverges.
`_tools/test_exponential_metric.py` holds the published $R_{rr}$ and the product $g_{tt}g_{rr}$ to those values.

## Step 8: the diagrams

Spacetime diagrams, at $m = 1$: the plane of $t$ and $r$ from $r = 0$ to $4m$ with the throat marked where $\partial_r(g_{\theta\theta}) = 0$; the line $y = z = 0$ of the Cartesian chart, with the throat at $x = \pm m$; the plane of $t$ and $R$ from the throat $R = e\,m$ out, the throat marked on its left edge; and the plane of $t$ and $u$.
`--verify` checks every ray against $ct \pm r_*$, with $r_*$ the tortoise coordinate above, which is the Curzon-Chazy particle's on its axis, where the plane's metric is the same.

Conformal diagrams: with $x = r_*(r) - r_*(m)$ over the whole line, $p, q = \arctan((ct \mp x)/\ell)$ at $\ell = 4m$ fill the diamond, the throat on its axis.
The right edges are null infinity and the left edges are $r = 0$, drawn as singular, on the strength of Step 7 and of Bronnikov, Fabris and Zhidenko.
The isotropic and harmonic charts cover the diamond and the areal chart its right half.

Embedding diagram: on a moment of $t$ the equator is $e^{2m/r}(dr^2 + r^2d\phi^2)$, the circle of $r$ has radius $\rho = re^{m/r}$, and $g_{rr} - (d\rho/dr)^2 = e^{2m/r}(m/r)(2 - m/r)$.
That is positive for $r > m/2$, where the surface stands in flat space, vertical on the throat and climbing as Flamm's paraboloid far out, and negative for $r < m/2$, where it is drawn in Minkowski space and tends to a light cone.
The two parts lie level on $r = m/2$, the circle of radius $e^2m/2$.
The harmonic, areal and Cartesian charts are checked to give the same circles and heights.
