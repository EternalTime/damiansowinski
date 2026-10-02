# Van Den Broeck's warp drive

Chris Van Den Broeck's drive of 1999 is Alcubierre's metric with the flat space of every slice multiplied by $B^2$.
This note records the three charts, what each is checked against, the pocket the diagrams declare, and why that pocket is Krasnikov's and not Van Den Broeck's own polynomial.

## Step 1: the Cartesian chart

Van Den Broeck's equation (1) is $ds^2 = -dt^2 + B^2(r_s)\left[(dx - v_s(t)f(r_s)\,dt)^2 + dy^2 + dz^2\right]$ with $c = 1$, $r_s = \sqrt{(x - x_s(t))^2 + y^2 + z^2}$ and $v_s = dx_s/dt$ (arXiv:gr-qc/9905084).
The chart publishes it with $c$ restored and with $f$ and $B$ left as arbitrary functions of the four coordinates, as Alcubierre's chart leaves $f$, so every component is an identity that assumes nothing about the profile.
`print_charts.van_den_broeck_alcubierre` checks that at $B = 1$ the metric is Alcubierre's published one, slot by slot.
The lapse is $1$ and the shift $(-v_sf, 0, 0)$, as in Alcubierre's, and the spatial metric is $B^2\delta_{ij}$, conformally flat and curved wherever $B$ varies.

## Step 2: the comoving spherical chart

Inside the wall of the bubble $f = 1$, and with $x' = x - x_s(t)$ the form $dx - v_s\,c\,dt$ is $dx'$ for any $v_s(t)$.
The metric there is $-c^2dt^2 + B^2(dx'^2 + dy^2 + dz^2)$, Van Den Broeck's $\mathrm{diag}(-1, B^2, B^2, B^2)$, static and spherically symmetric about the ship for $B = B(r)$, and the chart writes it in $(t, r, \theta, \phi)$ with $x' = r\cos\theta$.
`print_charts.van_den_broeck_pocket` pulls the Cartesian chart back through that map at $f = 1$, with $v_s = dx_s/d(ct)$ an arbitrary function, and compares every slot.
It also checks two of his equations against the chart's own tensors.
His energy density for the Eulerian observers, $T^{\hat 0\hat 0} = \frac{1}{8\pi}\left(B'^2/B^4 - 2B''/B^3 - 4B'/(B^3r)\right)$, is $G_{tt}/8\pi$.
His largest frame component of Riemann, $R_{\hat 1\hat 2\hat 1\hat 2} = B'^2/B^4 - B''/B^3 - B'/(B^3r)$ at a point of the line $y = z = 0$, where $e_1$ is radial and $e_2$ tangential, is $R_{r\theta r\theta}/(B^4r^2)$.

## Step 3: Krasnikov's chart

Krasnikov (arXiv:gr-qc/0207057, "Van Den Broeck's trick") writes the pocket as $ds^2 = -dt^2 + dl^2 + r(l)^2(d\theta^2 + \sin^2\theta\,d\phi^2)$ with $l \ge -l_0$ and $r(-l_0) = 0$, the Morris-Thorne line element with $\Phi = 0$ closed off at a centre.
With $\rho$ the comoving chart's radius, $dl = B\,d\rho$ and $r = B\rho$, which `print_charts.van_den_broeck_proper` checks by pulling the chart back.
The same function checks his Einstein tensor in the static frame, $G_{\hat t\hat t} = (1 - r'^2 - 2rr'')/r^2$, $G_{\hat t\hat t} + G_{\hat r\hat r} = -2r''/r$ and $G_{\hat t\hat t} + G_{\hat\theta\hat\theta} = (1 - r'^2 - rr'')/r^2$.
The weak energy condition fails exactly where $r'' > 0$, the part of the neck where the spheres shrink and then grow.

## Step 4: which pockets stand in flat space

On the equator of a slice the metric is $dl^2 + r(l)^2d\phi^2$, a surface of revolution in flat space with $dz/dl = \sqrt{1 - r'^2}$, which exists exactly where $|r'| \le 1$.
In the comoving chart $r' = d(B\rho)/(B\,d\rho) = 1 + \rho B'/B$, so the condition is $-2 \le \rho B'/B \le 0$.
Integrating $d\ln B/d\ln\rho \ge -2$ across the neck, from $\tilde R$ where $B = 1 + \alpha$ to $\tilde R + \tilde\Delta$ where $B = 1$, gives $1 + \alpha \le (1 + \tilde\Delta/\tilde R)^2$.
Van Den Broeck's example has $\tilde\Delta = \tilde R$ and $\alpha = 10^{17}$, far beyond the bound of $4$: across his neck the circles shrink faster than the distance in to them, and no surface of revolution in flat space carries the slice.
His polynomial $B = 1 + \alpha(n w^{n-1} - (n - 1)w^n)$ stays inside the bound for $\tilde\Delta = \tilde R$ only up to $\alpha$ of about $1$ at $n = 4$, where the neck is narrower than the pocket by less than a tenth.
Krasnikov's conditions on $r(l)$ include $|r'| \le 1$, so his pocket is the one that can be drawn, and the diagrams declare a pocket of his kind.

## Step 5: the declared pocket

In units of $R$, the radius where the wall of the bubble begins, and with $l = 0$ in the middle of the neck:

| part | $l$ | $r(l)$ | $r'$ |
| --- | --- | --- | --- |
| floor | $-4$ to $-5/2$ | $l + 4$ | $1$ |
| rim | $-5/2$ to $-3/2$ | $7/4 - (l + 2)^2$ | $1$ to $-1$ |
| lid | $-3/2$ to $-1/2$ | $-l$ | $-1$ |
| neck | $-1/2$ to $1/2$ | $l^2 + 1/4$ | $-1$ to $1$ |
| outside | from $1/2$ | $l$ | $1$ |

The neck is Krasnikov's own quadratic, $r = l^2/2l_1 + l_1/2$ with $l_1 = R/2$.
$r$ and $r'$ are continuous and $r''$ jumps at the four joins, so the curvature is bounded and discontinuous there, as it is for Krasnikov's.
Where $r' = \pm 1$ the slice is flat: the floor is a disc of radius $3R/2$, the lid the annulus from $3R/2$ in to $R/2$ described from outside in, and outside the neck the plane.
The surface climbs by $\int\sqrt{1 - 4u^2}\,du = \pi/4$ round the rim and by $\pi/4$ again through the neck, so the plane outside lies $\pi R/2$ above the floor.

The comoving radius follows from $d\ln\rho = dl/r$ with $\rho = l$ outside, and $B = r/\rho$:

| part | $\rho$ | $l(\rho)$ | $B(\rho)$ |
| --- | --- | --- | --- |
| outside | $\ge 1/2$ | $\rho$ | $1$ |
| neck | $\rho_c = e^{-\pi}/2$ to $1/2$ | $\tfrac12\tan(\tfrac12\ln 2\rho + \pi/4)$ | $1/(2\rho(1 - \sin\ln 2\rho))$ |
| lid | $\rho_c/3$ to $\rho_c$ | $-e^{-\pi}/4\rho$ | $e^{-\pi}/4\rho^2$ |
| rim | $\rho_0$ to $\rho_c/3$ | $-2 + \tfrac{\sqrt7}{2}\tanh\theta$ | $\tfrac{7}{4}\,\mathrm{sech}^2\theta/\rho$ |
| floor | $\le \rho_0$ | $-4 + \tfrac{3}{2}\rho/\rho_0$ | $3/2\rho_0$ |

with $\theta = \tfrac{\sqrt7}{2}\ln(3\rho/\rho_c) + \mathrm{artanh}(1/\sqrt7)$ and $\rho_0 = (\rho_c/3)\exp(-4\,\mathrm{artanh}(1/\sqrt7)/\sqrt7) = 0.003948$.
So $1 + \alpha = 3/2\rho_0 = 379.96$, and in Van Den Broeck's letters $\tilde R = 0.0039\,R$ and $\tilde R + \tilde\Delta = R/2$.
On the lid $B \propto \rho^{-2}$, the inversion $\rho \to 1/\rho$ of flat space.
`slices.checks` holds $B\rho = r(l)$, $B\,d\rho/dl = 1$ and the two closed forms of $l(\rho)$ and $\rho(l)$ to each other at four hundred points, and `embedding.van_den_broeck` checks that the comoving chart with this $B$ draws the same circles at the same heights and distances as Krasnikov's chart with this $r(l)$.

## Step 6: the diagrams

The wall of the declared bubble falls from $f = 1$ at $r_s = R$ to $0$ at $3R/2$ as $1 - 10s^3 + 15s^4 - 6s^5$ with $s = 2(r_s - R)/R$, twice differentiable, and $v_s = 2$ as Alcubierre's diagrams take it.
On the axis of the Cartesian chart the null condition is $dx/d(ct) = v_sf \pm 1/B$, and $\partial_y$ and $\partial_z$ of $f$ and $B$ vanish there, so the axis is totally geodesic and the curves drawn are null geodesics.
$g^{xx} = 1/B^2 - v_s^2f^2$ vanishes where $v_sfB = 1$; outside the neck that is Alcubierre's $v_sf = 1$.
In the comoving chart the radial rays are $ct \pm l(r)$ constant, and in Krasnikov's $ct \pm l$, both of which `null_rays.py --verify` checks.

No conformal diagram is drawn, for Alcubierre's reason: the causal structure of the whole spacetime is whatever $f$, $B$ and $v_s(t)$ make it.
The embedding diagram is the equator of a moment of $t$ in Krasnikov's chart, from the centre of the pocket to the circle $r = R$ where the wall begins.
Beyond it the slice stays flat through the wall, since the spatial metric is $B^2\delta_{ij}$ whatever $f$ is.
