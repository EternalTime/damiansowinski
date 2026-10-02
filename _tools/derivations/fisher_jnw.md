# The Fisher-Janis-Newman-Winicour scalar field

The static, spherically symmetric, asymptotically flat solution of Einstein's equations with a massless scalar field.
This note records the four charts, where each comes from, the field equations each is held to, and what each diagram draws.
Units here have $G = c = 1$ where a source does; the metric file keeps $c$ and writes the mass as a length.

## Step 1: the spherical chart

Virbhadra's (11) to (13): $ds^2 = -f^\gamma dt^2 + f^{-\gamma}dr^2 + f^{1-\gamma}r^2d\Omega^2$ with $f = 1 - b/r$, $\gamma = 2m/b$ and $b = 2\sqrt{m^2 + q^2}$, for a mass $m$ and a scalar charge $q$.
It is the form Agnese and La Camera gave Wyman's solution, and Claudel, Virbhadra and Ellis write it the same way in their (90) to (92), with $\nu$ for $\gamma$.
At $\gamma = 1$ it is Schwarzschild's metric with $r_s = b$.
For $\gamma < 1$ the sphere $r = b$ has the area $4\pi r^2f^{1-\gamma} = 0$ and the Kretschmann scalar diverges there as $(r - b)^{2\gamma - 4}$.
The metric file takes $b$ and $\gamma$ as the parameters, so the mass is $2GM/c^2 = \gamma b$ and the scalar charge, as a length, is $q = b\sqrt{1 - \gamma^2}/2$.

## Step 2: the field equations

With $R_{\mu\nu} = 2\,\partial_\mu\varphi\,\partial_\nu\varphi$ the field is $\varphi = \tfrac{1}{2}\sqrt{1 - \gamma^2}\,\ln f$, which is Virbhadra's $\Phi = (q/b\sqrt{4\pi})\ln f$ with his $R_{ij} = 8\pi\Phi_{,i}\Phi_{,j}$.
The only component of the Ricci tensor is $R_{rr} = b^2(1 - \gamma^2)/2r^2(r - b)^2$.
In the units of the metric file, $R_{\mu\nu} = (8\pi G/c^4)\,\partial_\mu\varphi\,\partial_\nu\varphi$ and $\varphi = \sqrt{(1 - \gamma^2)c^4/16\pi G}\,\ln f$.
`print_charts.fisher_jnw_check` holds every chart to $R_{\mu\nu} = 2\,\partial_\mu\varphi\,\partial_\nu\varphi$ in every slot, to the wave equation $\partial_x(\sqrt{-g}\,g^{xx}\partial_x\varphi) = 0$, written through logarithmic derivatives so that no root is taken, and to a vanishing Ricci tensor without the field.
The Kretschmann scalar is $K = b^2\left(48\gamma^2r^2 - 16\gamma(1 + \gamma)(1 + 2\gamma)br + (1 + \gamma)^2(3 + 2\gamma + 7\gamma^2)b^2\right)f^{2\gamma}/4r^4(r - b)^4$, which is $12b^2/r^6$ at $\gamma = 1$.

## Step 3: Janis, Newman and Winicour's radius, which is Fisher's function

Janis, Newman and Winicour's line element, Virbhadra's (5), has $g_{tt} = \left((1 - a_-/R)/(1 + a_+/R)\right)^{1/\mu}$ with $a_\pm = r_0(\mu \pm 1)/2$, and Virbhadra's (14) and (15) give $r = R + a_+$, $r_0 = 2m$ and $\mu = b/2m$.
So $R = r - b(1 + \gamma)/2$, $r_0 = \gamma b$, $\mu = 1/\gamma$, and $f = (2R - b(1 - \gamma))/(2R + b(1 + \gamma))$.
In this radius $g_{\theta\theta} = (R - a_-)^{1-\gamma}(R + a_+)^{1+\gamma}$, whose derivative is $2R\,g_{\theta\theta}/((R - a_-)(R + a_+))$, so $\Gamma^R{}_{\theta\theta} = -R$ for every $\gamma$.
Fisher wrote the metric in the areal radius, his (5), through $Z = re^{(\nu - \lambda)/2}$, his (14), and found $e^\nu = c^2\left((Z - Z_0)/(Z + Z_1)\right)^p$, his (25), with $Z_{0,1} = c^{-1}(\sqrt{(km)^2 + c^2a^2} \mp km)$ and $p = km/\sqrt{(km)^2 + c^2a^2}$, his (21).
Those are $Z_0 = c\,a_-$, $Z_1 = c\,a_+$ and $p = \gamma$, so $Z = cR$: Fisher's function is Janis, Newman and Winicour's coordinate.
The areal radius itself is given only implicitly, Fisher's (20), so it is no chart with components to print.

## Step 4: the isotropic radius

Fisher's (35) and (36) say that $Z = c\rho - km/c + (kG^2 + k^2m^2)/4c^3\rho$ brings the metric to isotropic form, which is $r = \rho(1 + b/4\rho)^2$.
Xanthopoulos and Zannias integrate the equations in that radius in any dimension, as Abdolrahimi and Shoom's note [21] records.
Then $f = h^2$ with $h = (4\rho - b)/(4\rho + b)$, and $ds^2 = -h^{2\gamma}dt^2 + (1 + b/4\rho)^4h^{2(1-\gamma)}(d\rho^2 + \rho^2d\Omega^2)$, with the singularity at $\rho = b/4$.
The reader holds $h^{2\gamma}$ as a power of $b - 4\rho$ with $(-1)^{2\gamma}$ beside it, a phase of its generator, so the chart's powers are counted on that base and the pullback is compared on positive symbols.

## Step 5: Bronnikov's harmonic coordinate

Bronnikov, Fabris and Zhidenko's (54), from Bronnikov's paper of 1973: $ds^2 = -e^{-2mu}dt^2 + \frac{k^2e^{2mu}}{\sinh^2(ku)}\left(\frac{k^2du^2}{\sinh^2(ku)} + d\Omega^2\right)$ with $\varphi \propto u$, and their (56), $e^{-2ku} = 1 - 2k/r$.
So $k = b/2$, $m = \gamma k$, spatial infinity is $u = 0$ and the singularity is $u \to \infty$; $u$ is an inverse length, which `DIMENSIONS` declares as `L**(-1)`.
$\sqrt{-g}\,g^{uu}$ does not depend on $u$, which is the harmonic condition and is checked, and $R_{uu} = 2(k^2 - m^2)$ is a constant.
With $k \to ik$ the same formulas give the wormholes of a field of the opposite sign, the Ellis-Bronnikov wormhole at $m = 0$.

## Step 6: the diagrams

Every diagram is drawn at $\gamma = 1/2$, the value of Abdolrahimi and Shoom's Figures 6, 9 and 11.
There $dr_*/dr = f^{-1/2}$ and $r_* = \sqrt{r(r - b)} + b\ln\left((\sqrt{r} + \sqrt{r - b})/\sqrt{b}\right)$, which vanishes at $r = b$; `null_rays.py --verify` checks the rays of all four charts against $ct \mp r_*$.
In the isotropic radius $d\rho/d(ct) = \pm h^{2\gamma - 1}(1 + b/4\rho)^{-2}$, which at $\gamma = 1/2$ is $\pm 1/4$ on the singularity, so that row stops its rays at the edge of the published domain.
The harmonic plane is drawn in units of $k$, at $k = 1$ and $m = 1/2$, where $b = 2k$.

Since $r_*$ is finite at $r = b$, $p, q = \arctan((ct \mp r_*)/\ell)$ bring the spacetime into Minkowski's triangle with the singularity a timelike line on $X = 0$, Abdolrahimi and Shoom's Figure 9; the drawing takes $\ell = 4b$, and one event is checked to land on one point through all four maps.

The equator of a moment of $t$ has circles of radius $\rho = rf^{(1-\gamma)/2}$ and $g_{rr} - (d\rho/dr)^2 = f^{-1-\gamma}\left(\gamma b/r - (1 + \gamma)^2b^2/4r^2\right)$, which vanishes on $r_e = (1 + \gamma)^2b/4\gamma$, where the Misner-Sharp energy $\tfrac{1}{2}\rho(1 - |\nabla\rho|^2)$ changes sign.
Outside $r_e$ the surface is drawn in flat space and inside it in three dimensional Minkowski space, Abdolrahimi and Shoom's (109) to (113).
At $\gamma = 1/2$, with $w = f^{1/4}$: $\rho = w/(1 - w^4)$, $r_e = 9b/8$ is $w^4 = 1/9$, and the height is the integral of $\sqrt{|9w^4 - 1|}/(1 - w^4)^{3/2}\,dw$ from $3^{-1/2}$, which `embedding.py` checks by a quadrature of its own.
The surface is read in the harmonic chart, where $w = e^{-u/4}$ and nothing cancels on the way to the singularity.
It enters the singular point along the light cone, with $1 - (dZ/d\rho)^2 = 16f$, so a chord's proper length is $4\sqrt{f}$ of its extent in $\rho$; the profile is written to fourteen decimals and stops at $u = 17.5$, a circle of radius $0.0126\,b$.
The geometry is static, so it is one surface and no movie.
