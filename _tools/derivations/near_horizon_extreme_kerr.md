# The near-horizon extreme Kerr geometry

James Bardeen and Gary Horowitz's throat of the extreme Kerr black hole, Phys. Rev. D 60, 104030 (1999), a vacuum with the isometry group $SL(2,\mathbb{R}) \times U(1)$ and one parameter, the angular momentum $J$.
This note records its four charts, the maps between them, and what each diagram draws.

## Step 1: the scale

Bardeen and Horowitz set $G = c = 1$ and write $r_0^2 = 2M^2$, with $J = M^2$ for the extreme hole.
With the constants restored the one length is $r_0$, $r_0^2 = 2GJ/c^3 = 2G^2M^2/c^4$, and the horizon's area $8\pi G^2M^2/c^4$ is $4\pi r_0^2$.
Guica, Hartman, Song and Strominger, Phys. Rev. D 80, 124008 (2009), write the same factor as $2GJ$ in front of the whole line element, with $\Omega^2 = (1 + \cos^2\theta)/2$ and $\Lambda = 2\sin\theta/(1 + \cos^2\theta)$.

## Step 2: the four charts

The Poincaré chart is Bardeen and Horowitz's throat metric as they write it, with $t$ a time and $r$ a length:
$ds^2 = \tfrac{1 + \cos^2\theta}{2}\left(-\tfrac{r^2}{r_0^2}c^2dt^2 + \tfrac{r_0^2}{r^2}dr^2 + r_0^2d\theta^2\right) + \tfrac{2r_0^2\sin^2\theta}{1 + \cos^2\theta}\left(d\phi + \tfrac{r}{r_0^2}c\,dt\right)^2$.
The inverse radius chart is the same patch in $x = r_0^2/r$, where the plane of $t$ and $x$ at fixed angles is conformal to flat: Guica and her coauthors write it with the pure number $y = x/r_0$, and Bardeen and Horowitz use $x = 1/r$ to bring the boundary to a finite distance.
The global chart is Bardeen and Horowitz's, with $\tau$ and $y$ pure numbers and $r_0^2$ in front: $-(1 + y^2)d\tau^2 + dy^2/(1 + y^2)$ in place of the Poincaré bracket and $d\phi + y\,d\tau$ in the twist.
The near-NHEK chart is the throat in coordinates that end at a horizon of finite temperature, Amsel, Horowitz, Marolf and Roberts, JHEP 09 (2009) 044, with $r(r - 2k)/r_0^2$ in place of $r^2/r_0^2$ and $(r - k)/r_0^2$ in the twist; Bredberg, Hartman, Song and Strominger, JHEP 04 (2010) 019, named it near-NHEK and put its horizon at $r = 0$, which is this chart with $r$ moved by $2k$.
Its surface gravity with respect to $ct$ is $f'(2k)/2 = k/r_0^2$ for $f = r(r - 2k)/r_0^2$, which is what the parameter's description states.

## Step 3: the maps

`print_charts.near_horizon_extreme_kerr_pullback` pulls the Poincaré chart back through each map and compares it with the chart's own metric in every slot.
The inverse radius is $r = r_0^2/x$ with the same $t$ and $\phi$.
The global chart enters by Bardeen and Horowitz's map with $r_0$ restored: $r = r_0\left(\sqrt{1 + y^2}\cos\tau + y\right)$, $ct = r_0\sqrt{1 + y^2}\sin\tau/\left(\sqrt{1 + y^2}\cos\tau + y\right)$ and $\phi \to \phi + \ln\left((\cos\tau + y\sin\tau)/(1 + \sqrt{1 + y^2}\sin\tau)\right)$.
On $\tau = 0$ that is $t = 0$, $r = r_0\left(\sqrt{1 + y^2} + y\right)$ and the same $\phi$, which is how the moment the embedding diagram draws is carried onto the Poincaré planes.
With $y = \tan\sigma$ the Poincaré null coordinates are $ct/r_0 + x/r_0 = \tan\left((\tau - \sigma + \pi/2)/2\right)$ and $ct/r_0 - x/r_0 = \tan\left((\tau + \sigma - \pi/2)/2\right)$, so the Poincaré wedge sits in the strip of $\tau$ and $\sigma$ as `conformal.poincare_pq` puts it.
The near-NHEK chart enters by the exponentials of its null coordinates: with $r_* = (r_0^2/2k)\ln(1 - 2k/r)$, $U = e^{k(ct + r_*)/r_0^2}$ and $V = e^{k(ct - r_*)/r_0^2}$ are the Poincaré chart's $ct/r_0 - r_0/r$ and $ct/r_0 + r_0/r$, and $\phi$ moves by $\tfrac{1}{2}\ln(1 - 2k/r)$, since $(r - k)\,c\,dt/r_0^2 - r_P\,c\,dt_P/r_0^2 = (k/r_0^2)\,dr_*$ is exact.
The patch is $U > 0$, the part of the Poincaré wedge to the future of the ray $ct = r_0^2/r$, and its future horizon $V \to \infty$ is a stretch of the Poincaré horizon.

## Step 4: the curvature

All four charts are Ricci flat, which the script checks before it writes, so the Weyl tensor is the Riemann tensor.
The Kretschmann scalar is Kerr's, $48M^2(r^2 - a^2\cos^2\theta)\left((r^2 + a^2\cos^2\theta)^2 - 16a^2r^2\cos^2\theta\right)/(r^2 + a^2\cos^2\theta)^6$, at $r = a = M$: $K = 192\left(\cos^4\theta - 14\cos^2\theta + 1\right)\sin^2\theta/r_0^4(1 + \cos^2\theta)^6$, a function of $\theta$ alone, $192/r_0^4$ on the equator and zero on the axis.
Every value is printed in $1 + \cos^2\theta$, the factor the line element is written in, and the near-NHEK chart's in powers of $r - k$ and $k$.

## Step 5: the rays

The two null directions of no angular momentum, $k = (1, \pm v, 0, w)$ with $w = -g_{0\phi}/g_{\phi\phi}$, are geodesic and are repeated principal null directions of the Weyl tensor at every $\theta$, which `print_charts.near_horizon_extreme_kerr_checks` checks in each chart.
Their shadows on the plane of the time and the radius are the null curves of the metric orthogonal to the circles of $\phi$, $\Omega^2$ times a two dimensional anti-de Sitter space of radius $r_0$, the same curves at every $\theta$.
So the spacetime diagrams draw the equator with $\phi$ divided out, `quotient="phi"`, and the conformal diagram is the strip of that anti-de Sitter space.
On the equator $g_{tt}$ of the Poincaré chart is $3r^2/2r_0^2$, positive at every $r$, and $g_{\tau\tau}$ of the global chart is $r_0^2(3y^2 - 1)/2$, which vanishes at $y = \pm 1/\sqrt{3}$.
Off the equator $\partial_t$ of the Poincaré chart is spacelike where $\cos^4\theta + 6\cos^2\theta < 3$, which is $\cos^2\theta < 2\sqrt{3} - 3$, within $42.9°$ of the equator; Bardeen and Horowitz's text gives $32.4°$, which is the angle whose sine is their $0.536$, and the condition they write, $\sin\theta > (1 + \cos^2\theta)/2$, has the root $\sin\theta = \sqrt{3} - 1 = 0.732$.

## Step 6: the embedding

The equator at the moment $\tau = 0$ is $r_0^2\,dy^2/2(1 + y^2) + 2r_0^2\,d\phi^2$, a cylinder of radius $\sqrt{2}\,r_0 = 2GM/c^2$ along which $z = r_0\,\mathrm{arsinh}(y)/\sqrt{2}$.
The sphere of $\theta$ and $\phi$ at any event is the horizon of the extreme hole, $g_{\theta\theta} = r_0^2(1 + \cos^2\theta)/2$ and $\rho = \sqrt{2}\,r_0\sin\theta/\sqrt{1 + \cos^2\theta}$.
$g_{\theta\theta} - (d\rho/d\theta)^2 = r_0^2\left((1 + c^2)^4 - 16c^2\right)/2(1 + c^2)^3$ with $c = \cos\theta$, and $(1 + c^2)^4 - 16c^2 = (c^2 - 1)(c^3 - c^2 + 3c + 1)(c^3 + c^2 + 3c - 1)$, negative from each pole to the root $c_1 = 0.2956$ of $c^3 + c^2 + 3c = 1$, $\theta_1 = 72.81°$, and positive between.
So the belt within $17.19°$ of the equator stands in flat space and each cap in three dimensional Minkowski space, where it climbs at $dZ/d\theta = \sqrt{(d\rho/d\theta)^2 - g_{\theta\theta}}$; on $\theta_1$ the surface lies level in both, and `embedding.LevelSlice` states that tangent exactly, since the root has no form sympy reduces.
The Gaussian curvature is $4(1 - 3\cos^2\theta)/r_0^2(1 + \cos^2\theta)^3$, negative within $54.74°$ of each pole, the polar caps Smarr found on a Kerr horizon with $a > \sqrt{3}\,GM/2c^2$, Phys. Rev. D 7, 289 (1973).
