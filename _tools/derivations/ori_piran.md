# Self similar fluid collapse to a naked singularity

## Sources and charts

Both charts are Tomohiro Harada and Hideki Maeda's, arXiv:gr-qc/0101064, section II, in geometrized units $G = c = 1$.
The comoving chart is their line element (1) with every dimensionless function of $z = r/(-t)$ alone, as their section II.B assumes; Carr and Gundlach, arXiv:gr-qc/0209092, section III, write the same form with $z = r/t$ and the time reversed.
The closed forms $e^{\sigma} = a_\sigma z^{4k/(1+k)}\eta^{-2k/(1+k)}$ and $e^{\omega} = a_\omega\eta^{-2/(1+k)}S^{-4}$ follow their equations (dsigmadr) and (domegadt), with $a_\omega = 1$ and the gauge $a_\sigma = (2D)^{2k/(1+k)}$ that sets $e^{\sigma} = 1$ at the centre, so $e^{\sigma} = (2Dz^2/\eta)^{2k/(1+k)}$.
The zooming chart is their $\tau = -\ln(-t)$ with a length $\ell$, $t = -\ell e^{-\tau}$ and $r = \ell ze^{-\tau}$, which puts the scaling symmetry in $\partial_\tau$; its check holds it to the comoving metric pulled back.
In the comoving chart $\sigma$, $\omega$ and $S$ are written as functions of $t$ and $r$, since the checker reads a function of the chart's own coordinates, and their descriptions say they depend on $z$ alone; the zooming chart writes them as functions of its coordinate $z$.

## The equations and the boundary conditions

In $x = \ln z$, with $y = M/(\eta S^3)$ and $V^2 = z^{2(1-k)/(1+k)}\eta^{-2(1-k)/(1+k)}S^{-4}/a_\sigma$, Harada and Maeda's (dmdz), (dsdz) and their dependent equation for $\eta$ are

$$(\ln M)' = \frac{k}{1+k}\left(\frac1y - 1\right),\qquad (\ln S)' = -\frac{1-y}{1+k},\qquad (\ln\eta)' = 2 - \frac{4yV^2 - (1+k)^2S^{-4}\eta^{-(1-k)/(1+k)}}{2(V^2 - k)},$$

under the constraint (constraint), $(1-y)^2V^2 - (k+y)^2 + (1+k)^2\eta^{-2/(1+k)}S^{-6}(1 - M/S) = 0$.
The comoving chart's check holds its published Einstein tensor to $8\pi$ times the stresses of the fluid at rest, $\rho = \eta/(8\pi r^2)$ and $P = k\rho$, at two states on the constraint surface with $\ln S$ and $\ln\eta$ carried to second order along these equations, which tests the equations as transcribed against Einstein's.
`_tools/test_ori_piran.py` does the same with the published tensor on the numerical solution itself, at five values of $z$, to $2\times10^{-5}$.

A regular centre has Harada and Maeda's series, $y \to 1/3$ and $\eta \simeq 2Dz^2$, with $D = 4\pi\rho_c t^2$.
The sonic point is $V^2 = k$; there the numerator of $(\ln\eta)'$ must vanish as well, and the two values of $(\ln\eta)'$ through a point of the sonic line are the roots of a quadratic from l'Hôpital's rule, `sonic_slopes`.

## The numerical solution

`ori_piran.py` integrates the relativistic Larson-Penston solution at $k = 1/100$.
The centre's series is integrated outward from $z = 10^{-12}$, and the analytic expansion at a point of the sonic line, along the smaller of its two slopes, inward; Newton's method on $D$ and on $y$ at the sonic point makes the two agree in $\ln S$ and $\ln\eta$ at $z = 10^{-3}$, and the same expansion is integrated outward to $z = 10^8$.
Three starting guesses reach the same solution, and steps of $10^{-5}$ to $10^{-7}$ off the sonic line and matching points of $10^{-3}$ and $5\times10^{-3}$ change $D$ in its eighth figure.
The result is $D = 1.4377137$, $y = 0.4145670$ at the sonic point, the sonic point at $z_s = 0.0099012$ and the null orbit of the scaling symmetry, $V_z^2 = 1$, at $z_1 = 0.2013946$.
Harada and Maeda give $D \approx 1.439$ at $k = 0.01$ and say that it coincides with Ori and Piran's; the integration here differs from it by $0.09$ per cent.
These are computed diagnostics, not numbers taken from a paper; the History prints $D \approx 1.44$, the value both round to, with Harada and Maeda's citation.
Matching along the larger slope did not converge from the same starting guesses.

Along the solution $2m/R < 1$ and $e^{-\omega}(\partial_rR)^2 < 1$ for every $t < 0$, so no sphere is trapped before the singularity forms and the slice of constant $t$ embeds in flat space.
Far out $y \to 1$ and both rates fall as $z^{-q}$, with $q = 0.98015$ from the last grid points, which is $(1-k)/(1+k)$; past the grid the tables continue with that decay, and below it with the centre's limiting rates.
The tables hold $\ln S$, $\ln\eta$ and their exact rates from the equations, each a cubic spline in $x$, so every first derivative is a rate of the equations and every second the derivative of one; differentiating the logarithms instead loses the rates to rounding once they fall below $10^{-14}$.

## Drawings

The spacetime diagrams draw the comoving plane of $t$ and $r$ over $0 \le r \le 0.7\ell$, $-1.2\ell \le t \le 0$, and the zooming plane of $\tau$ and $z$ over $0 \le z \le 0.3$, $-0.1 \le \tau \le 0.2$, each holding the moment $t = -\ell$ of the embedding diagram inside it.
Both mark the sonic point and the null orbit, lines of constant $z$ through the scaling origin: by Carr and Gundlach's section III a line of constant $z$ with $V = \pm1$ is a radial null geodesic, and the one with $t < 0$ that reaches $r = 0$ at $t = 0$ is the past light cone of the singularity; Guo, Hadzic and Jang's sonic line is the boundary of the past sound cone of the scaling origin.
Neither marks where the Kretschmann scalar diverges: the only such point in either window is the corner $t = r = 0$, and the published scalar, expanded, cancels to no digits in doubles where $z$ is large, so it would mark the regular moment $t = 0$ falsely.

The conformal diagram is not drawn.
Both sourced charts end at $t = 0$, where the singularity forms, and the light that escapes it runs in the continuation past $t = 0$, which Guo, Hadzic and Jang construct by a maximal analytic extension in their adapted comoving coordinates and which no chart here covers.
Their adapted comoving chart is not given: before $\tilde\tau = 0$ its functions are the comoving ones rescaled, their (E:MUTILDE) and (E:LAMBDATILDE), so it would repeat the comoving chart, and its continuation needs that extension, which the integration here does not carry.

The embedding diagram is the comoving slice $t = -\ell$ from $r = 10^{-6}\ell$ to $0.4\ell$, where $e^{-\omega}(\partial_rR)^2$ falls from 1 to 0.940, with rings at $z_s$ and $z_1$.
The curvature figure plots the stored $K$ and $R$ along the same slice from $r = 10^{-3}\ell$, nearer to the centre than which the expanded $K$ disagrees with its sixty-digit value by more than $10^{-9}$, to $r = \ell$; it checks $R = 8\pi(1 - 3k)\rho$ against the solution's own $\eta$, the trace of the field equations, and that $8\pi(1-3k)\rho$ tends to $2(1-3k)D$ at the centre.
