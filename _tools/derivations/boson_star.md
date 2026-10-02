# Boson stars and geons

A boson star is a static, spherically symmetric ball of complex scalar field held together by its own gravity.
No closed form is known for its metric.
The entry therefore gives two charts with their metric functions left free, as `tov` does, and every drawing is made from one numerical solution, which `_tools/derivations/boson_star.py` computes and `_tools/test_boson_star.py` holds to the physics.

The charts are written by `_tools/derivations/print_charts.py --metric boson_star`, and `verify_metrics.py --system boson_star/<chart>` checks each in about a second.

## Step 1: the charts and their sources

- `areal`: S. L. Liebling and C. Palenzuela, Living Rev. Relativ. 26, 1 (2023), arXiv:1202.5809v5, their (33), "polar-areal coordinates",

  $$ds^2 = -\alpha(r)^2c^2dt^2 + a(r)^2dr^2 + r^2d\Omega^2.$$

  It is the chart their equations (37) to (39) are integrated in, and C.-W. Lai's thesis, arXiv:gr-qc/0410040, has the same line element as his (4.90).
  The mass inside $r$ is $M(r) = c^2r(1 - a^{-2})/2G$, their (46).
- `isotropic`: Liebling and Palenzuela's (47),

  $$ds^2 = -\alpha(R)^2c^2dt^2 + \psi(R)^4\left(dR^2 + R^2d\Omega^2\right),$$

  reached from the first by $dR/dr = aR/r$ and $\psi = \sqrt{r/R}$, their (49) and Lai's Appendix D.
  F. E. Schunck and E. W. Mielke, Class. Quantum Grav. 20, R301 (2003), arXiv:0801.0307, print their line element (2.7) in the same isotropic form, with $e^{\nu}$ and $e^{\lambda}$ for $\alpha^2$ and $\psi^4$.

Kaup's and Ruffini and Bonazzola's own papers were reached only through their abstracts, so the notation is the review's.
Wheeler's idealized spherical geon has a static, spherically symmetric metric as well, the time average of its field's, and the areal chart with its functions free covers it; no geon is drawn.

## Step 2: the field equations

The action is Liebling and Palenzuela's (7) and (8), $\mathcal{L}_M = -\tfrac{1}{2}\left[g^{ab}\nabla_a\bar\varphi\nabla_b\varphi + \mu^2|\varphi|^2\right]$, with the stress tensor of their (10).
For $\varphi = \varphi_0(r)e^{i\omega t}$ and $\sigma = \sqrt{4\pi G}\,\varphi_0/c^2$, in units with $c = 1$,

$$\frac{8\pi G}{c^4}T^t{}_t = -\left(\frac{\omega^2\sigma^2}{\alpha^2} + \frac{\sigma'^2}{a^2} + \mu^2\sigma^2\right),\qquad \frac{8\pi G}{c^4}T^r{}_r = \frac{\omega^2\sigma^2}{\alpha^2} + \frac{\sigma'^2}{a^2} - \mu^2\sigma^2,\qquad \frac{8\pi G}{c^4}T^\theta{}_\theta = \frac{\omega^2\sigma^2}{\alpha^2} - \frac{\sigma'^2}{a^2} - \mu^2\sigma^2.$$

The stress along the radius differs from the stress across it by $2\sigma'^2/a^2$, which is why the star has no equation of state.
The sign of $\omega$ does not enter, so Lai's and Schunck and Mielke's $e^{-i\omega t}$ gives the same star.

With the published $G^t{}_t = -(a^3 - a + 2r\,\partial_r a)/(a^3r^2)$ and $G^r{}_r = (2r\,\partial_r\alpha + \alpha - \alpha a^2)/(\alpha a^2r^2)$ the first two give Liebling and Palenzuela's (37) and (38), and the wave equation gives their (39):

$$\partial_r a = \frac{a}{2}\left[-\frac{a^2 - 1}{r} + r\left(\left(\frac{\omega^2}{\alpha^2} + \mu^2\right)a^2\sigma^2 + \sigma'^2\right)\right],$$

$$\partial_r\alpha = \frac{\alpha}{2}\left[\frac{a^2 - 1}{r} + r\left(\left(\frac{\omega^2}{\alpha^2} - \mu^2\right)a^2\sigma^2 + \sigma'^2\right)\right],$$

$$\sigma'' = -\left(1 + a^2 - \mu^2r^2a^2\sigma^2\right)\frac{\sigma'}{r} - \left(\frac{\omega^2}{\alpha^2} - \mu^2\right)a^2\sigma.$$

`boson_star_check` in `print_charts.py` holds the chart to them in sympy: with $\partial_r a$ and $\partial_r\alpha$ solved from the published $G^t{}_t$ and $G^r{}_r$ and $\sigma''$ from the wave equation, the published $G^\theta{}_\theta$ less the field's stress across the radius is zero identically, which is the contracted Bianchi identity.
No Einstein component stands off the diagonal, and at $\sigma = 0$ with $a^{-2} = \alpha^2 = 1 - r_s/r$ the chart is the published metric of `schwarzschild`.
The isotropic chart is held to being the first pulled back along $r = \psi^2R$ with $a = \psi/(\psi + 2R\,\partial_R\psi)$.

## Step 3: the numerical solution

`boson_star.Star` works in units with $G = c = \mu = 1$, so a length is counted in $\hbar/mc$ and a mass in $M_P^2/m$.
It reads the published $G^t{}_t$ and $G^r{}_r$ from the metric file, as `null_rays.StarSolver` does for `tov`, and leaves the published $G^\theta{}_\theta$ for a check.

A solution regular at the centre has $a(0) = 1$, $\sigma'(0) = 0$ and a chosen $\sigma(0)$, and falls off far away only for a discrete set of frequencies.
The lowest, the state with no node, is found by bisection on $\omega^2$ with $\alpha(0) = 1$: too low and the field turns back up before it reaches zero, too high and it crosses zero.
This is the shooting of Liebling and Palenzuela's section 2.4 and of Lai's section 4.3.
The field of the last frequency that turns up is followed to its turning point, the solver's `edge`, at $\mu r = 39.7$ for the declared star, where $\sigma = 6\times10^{-10}$.
From there on the metric is Schwarzschild's of the mass inside, which differs from the star's by the square of what is left of the field, and $\alpha$ is scaled to meet it, which scales $\omega$ with it, their (44) and (45).

Within $0.1/\mu$ of the centre the three functions are their power series in $r^2$, `centre_series`, and the integration starts from the series there.
The equations divide differences that vanish as $r^2$ by $r$, so an integration started at the centre loses the second derivatives of the metric there, and with them the curvature; with the series $(a - 1)/r^2$ and $a'/r$ are finite sums.

The isotropic radius is $R = r\exp(L(r))$ with $L' = (a - 1)/r$, fixed at the edge by Schwarzschild's $R = (r - M + \sqrt{r^2 - 2Mr})/2$.
Liebling and Palenzuela's (48) and Lai's Appendix D write that value as $((1 + \sqrt{a})/2)^2r/a$, which holds when $a$ stands for $g_{rr}$; in the metric function $a$ of (33) it is $((1 + a)/2)^2r/a^2$.

## Step 4: the declared star

Every drawing declares the ground state with $\sigma(0) = 0.271$, the heaviest:

| | |
| --- | --- |
| mass $M$ | $0.633\,M_P^2/m$ |
| frequency $\omega$ | $0.853\,\mu c$ |
| lapse at the centre | $0.686$ |
| number of bosons $N$ | $0.653\,M_P^2/m^2$ |
| radius holding 99% of the mass, $R_{99}$ | $7.86/\mu$ |
| largest $2GM(r)/c^2r$ | $0.240$, at $\mu r = 3.82$ |
| conformal factor at the centre | $1.168$ |

$0.633\,M_P^2/m$ is Kaup's limit as Liebling and Palenzuela's section 2.4 and Schunck and Mielke's (2.19) quote it, Schunck and Mielke's Table II has the particle number $0.653$, and C. A. R. Herdeiro, A. M. Pombo and E. Radu, arXiv:1708.05674, Table 1, have $\omega = 0.853$ at the largest mass.

`_tools/test_boson_star.py` holds the solver to:

- the published $G^\theta{}_\theta$, to $10^{-10}$ of the central energy density;
- Kaup's limit, the largest mass of the family at the declared star;
- the first law $dM = \omega\,dN$ along the family, and $M < mN$ below the limit;
- the Newtonian limit, where a quarter of the central field gives half the mass, a quarter of $1 - \omega$ and twice the radius;
- a regular centre, where the series and the integration agree to ten places;
- Schwarzschild's metric from the edge on;
- the isotropic chart being the same star, with $r = \psi^2R$, $a = \psi/(\psi + 2R\psi')$ and the mass read off $\psi \to 1 + M/2R$.

It holds the drawings to the solver from the numbers in their files: every ray of both charts is a null curve, $c\,dt = \pm(a/\alpha)\,dr$, the marked sphere holds 99% of the mass, the embedded surface climbs at $\sqrt{a^2 - 1}$, and the conformal diagram lies in Minkowski's triangle.

## Step 5: dimensions

`DIMENSIONS` declares $t$ a time, $r$ and $R$ lengths, and $\alpha$, $a$ and $\psi$ pure numbers.
With $c$ and $G$ restored the frequency enters as $\omega^2/c^2\alpha^2$ beside $\mu^2$, both inverse areas, and $\sigma$ is dimensionless, as the parameters' descriptions write the three equations.
