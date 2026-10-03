# Flat space with supertranslation hair

Geoffrey Compère and Jiang Long, "Vacua of the gravitational field", JHEP 2016(07), 137, arXiv:1601.04958, and "Classical static final state of collapse with supertranslation memory", Class. Quantum Grav. 33, 195001 (2016), arXiv:1602.05197.
Minkowski space carried by a finite supertranslation with the field $C$ on the sphere, in their static coordinates (1601.04958, section 2.2)

$$ds^2 = -dt^2 + d\rho^2 + \left[(\rho - C)^2\gamma_{AB} - 2(\rho - C)D_AD_BC + D_AD_EC\,D_BD^EC\right]dz^Adz^B,$$

with the constant part of $C$ left out, and in BMS gauge about future null infinity (section 2.1) and about past null infinity (section 2.2).
Its three charts are written by `_tools/derivations/print_charts.py --metric supertranslation_hair`, and `verify_metrics.py --system supertranslation_hair/<chart>` checks each in seconds.

## Step 1. The field taken axisymmetric

Every chart takes $C = C(\theta)$, a length.
The Hessian of $C$ on the unit sphere is then diagonal, $D_\theta D_\theta C = C''$ and $D_\phi D_\phi C = \sin\theta\cos\theta\,C'$, so the metric of a shell is $(\rho - C - C'')^2d\theta^2 + B^2d\phi^2$ with $B = (\rho - C)\sin\theta - C'\cos\theta$.
Both brackets are $M_{AB} = (\rho - C)\gamma_{AB} - D_AD_BC$ squared through $\gamma^{AB}$, which is how the general line element factors.
A constant in $C$ moves the zero of $\rho$, and $a\cos\theta$ moves the origin by $a$ along the axis (1602.05197, section 3.1).

## Step 2. The charts and their sources

- `static`: their static form, 1601.04958 section 2.2, with $C(\theta)$, which `supertranslation_hair_check` holds to the flat metric pulled back through the Cartesian coordinates of 1602.05197 section 3.1, $X = (\rho - C)\,n - C'\,e_\theta$ for axisymmetric $C$, $n$ the unit vector of the angles and $e_\theta = \partial_\theta n$.
- `bondi_retarded`: their BMS gauge, 1601.04958 section 2.1 with the radius of section 2.2, $u = t - \rho$ and $\rho = \sqrt{r^2 + U} + \tfrac{1}{2}(D^2 + 2)C$ with $C_{AB} = -(2D_AD_B - \gamma_{AB}D^2)C$ and $U = C_{AB}C^{AB}/8$. For $C(\theta)$, $C_{\theta\theta} = 2\sigma$ and $C_{\phi\phi} = -2\sigma\sin^2\theta$ with the shear $\sigma = \tfrac{1}{2}(C'\cot\theta - C'')$, so $U = \sigma^2$, $\sqrt{r^2 + U} = W$, and the angular metric $(r^2 + 2U)\gamma_{AB} + W C_{AB}$ is $(W + \sigma)^2d\theta^2 + (W - \sigma)^2\sin^2\theta\,d\phi^2$. With $P = C + \tfrac{1}{2}(C'' + C'\cot\theta)$ one has $\partial_\theta P = -\partial_\theta\sigma - 2\sigma\cot\theta$, so $g_{u\theta} = -\partial_\theta(W + P)$ holds $C$ through $\sigma$ alone. `supertranslation_hair_check` holds the chart to the static chart pulled back at two fields $C$, and to Bondi's axisymmetric metric as published in `bondi_sachs.json` at $e^{2\beta} = r/W$, $re^\gamma = W + \sigma$, $U_{\rm Bondi} = -g_{u\theta}/(W + \sigma)^2$ and $V = W(1 + (W + \sigma)^2U_{\rm Bondi}^2)$, so $\sigma$ is the shear of that page.
- `bondi_advanced`: their advanced coordinates, 1601.04958 section 2.2, $v = t + \rho = u + 2\rho$, the retarded chart with $u \to -v$, held to the static chart pulled back.

## Step 3. The held names

$B$ in the static chart and $W$ in the Bondi charts are held, with their rates in `RATES` of `verify_metrics.py`: $\partial_\rho B = \sin\theta$, $\partial_\theta B = (\rho - C - C'')\cos\theta$, $\partial_rW = r/W$ and $\partial_\theta W = \sigma\,\partial_\theta\sigma/W$.
The rates fix each name up to a constant, and the metric is flat for every value of it, which was checked numerically for $W^2 = r^2 + \sigma^2 + k$, so no value holds a relation the rates leave out.
`Reader.by_rates` leaves the derivatives of $C$ and $\sigma$, declared functions and not held names, as they stand.

## Step 4. The supertranslation horizon

The determinant of the shells' metric vanishes where $\rho = C + \max(C'', C'\cot\theta)$, $\rho_{SH} = E + \sqrt{U}$ of 1602.05197 section 3.2; in the Bondi charts that is $r = 0$.
The drawings take $C = \tfrac{1}{2}\ell(3\cos^2\theta - 1)$ at $\ell = 1$, so $\sigma = -\tfrac{3}{2}\ell\sin^2\theta$; on the equator the horizon is $\rho = 5\ell/2$, the ring of radius $3\ell$ in flat space, and there $W = 3\ell/2$ at $r = 0$, $P = \ell$, so $t = u + W + \ell = v - W - \ell$ and $\rho = W + \ell$.
Each shell is convex, with principal radii $\rho - \ell(5 - 9\cos^2\theta)/2$ and $B/\sin\theta = \rho + \ell(1 + 3\cos^2\theta)/2$, $2(\rho + \ell/2)$ across the equator and $2(\rho - \ell)$ from pole to pole.
`_tools/test_supertranslation_hair.py` holds these numbers.
