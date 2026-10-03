# Kerr black holes with scalar hair

Carlos Herdeiro and Eugen Radu's rotating black holes with a cloud of massive complex scalar field are solutions of the Einstein-Klein-Gordon equations known only numerically.
The entry therefore gives one chart with its four metric functions left free, as `neugebauer_meinel` does, and the member with no field, Kerr's metric in the same coordinates, in closed form.
Every drawing is made from one solution, the authors' configuration IV, which `_tools/derivations/kerr_scalar_hair.py` reads from the data they published and `_tools/test_kerr_scalar_hair.py` holds to the physics.

The charts are written by `_tools/derivations/print_charts.py --metric kerr_scalar_hair`, and `verify_metrics.py --system kerr_scalar_hair/<chart>` checks each.

## Step 1: the charts and their sources

- `herdeiro_radu`: C. A. R. Herdeiro and E. Radu, Phys. Rev. Lett. 112, 221101 (2014), arXiv:1403.2757, their (3), and C. Herdeiro and E. Radu, Class. Quantum Grav. 32, 144001 (2015), arXiv:1501.04319, their (2.5),

  $$ds^2 = e^{2F_1}\left(\frac{dr^2}{N} + r^2d\theta^2\right) + e^{2F_2}r^2\sin^2\theta\left(d\varphi - W\,c\,dt\right)^2 - e^{2F_0}N\,c^2dt^2,\qquad N = 1 - \frac{r_H}{r},$$

  with $F_0$, $F_1$, $F_2$ and $W$ functions of $r$ and $\theta$.
  The paper writes $W\,dt$ in units with $c = 1$; with $x^0 = ct$ the chart's $W$ is an inverse length, and the horizon turns at $\Omega_H = cW(r_H)$.
  At $r_H = 0$ it is the ansatz of spinning boson stars, the paper's section 4.1.
- `kerr_member`: the paper's appendix A, (A.1), Kerr's metric in the same coordinates with the constants $r_H$ and $c_t < 0$.
  `_tools/derivations/kerr_scalar_hair.py` and the chart check confirm that it is Boyer and Lindquist's chart with $R = r - c_t$, $GM/c^2 = r_H/2 - c_t$ and $a^2 = c_t(c_t - r_H)$.
  So $-c_t$ is the Boyer-Lindquist radius of the inner horizon and $r_H = R_+ - R_-$, which the paper's $r = R - a^2/R_H$ says too, since $a^2/R_+ = R_-$.
  The chart is written in $b = -c_t > 0$: the checker's canonical form splits a root into the roots of its factors as if each were positive, and $\sqrt{c_t(c_t - r_H)}$ has two negative factors, while $\sqrt{b(b + r_H)}$ has none.

## Step 2: the field equations

The field is $\Psi = \phi(r,\theta)e^{i(m\varphi - wt)}$, the paper's (2.6), and the stress tensor is the paper's (2.2).
With the amplitude scaled as in the paper's (3.21), the $8\pi G\phi^2$ of its equations is $2\sigma^2$ for $\sigma = \sqrt{4\pi G}\,\phi/c^2$, the dimensionless amplitude `boson_star` uses.

`kerr_scalar_hair_check` in `print_charts.py` builds the field's stress from the chart's own metric and holds the chart's Einstein tensor to the paper's equations exactly:

- the wave operator of the chart, $(\Box - \mu^2)\Psi$, is the paper's (2.10) times $e^{-2F_1}N$;
- $E^r{}_r + E^\theta{}_\theta - E^\varphi{}_\varphi - E^t{}_t$, $E^r{}_r + E^\theta{}_\theta - E^\varphi{}_\varphi + E^t{}_t + 2WE^t{}_\varphi$, $E^r{}_r + E^\theta{}_\theta + E^\varphi{}_\varphi - E^t{}_t - 2WE^t{}_\varphi$ and $E^t{}_\varphi$ of $E^a{}_b = G^a{}_b - 2S^a{}_b$ are the combinations of its (2.11), and give the equations for $F_1$, $F_2$, $F_0$ and $W$ it prints after it, each times one factor;
- $E^r{}_r - E^\theta{}_\theta$ and $E^r{}_\theta$ are its constraints (2.12) and (2.13).

The paper's $E_\varphi{}^t$ has its $t$ up, $E^t{}_\varphi$; read with the $t$ down the combinations mix the equation for $W$ into those for $F_2$ and $F_0$.

## Step 3: the solution

`kerr_scalar_hair_IV.dat.gz` is the authors' configuration-IV.dat from http://gravitation.web.ua.pt/node/416, the paper's reference [gravwebsite], gzipped and otherwise unchanged; the test holds its SHA-256.
It tabulates $F_1$, $F_2$, $F_0$, $\phi$ and $W$ on 251 points of $X = x/(1 + x)$, $x = \sqrt{r^2 - r_H^2}$, from 0 to 1, by 30 of $\theta$ from 0 to $\pi/2$, in units $G = c = \mu = 1$, for $m = 1$, $w = 0.82$ and $r_H = 0.1$.
`Hair` reads it into bicubic splines in $X$ and $\theta$, with the angle continued through the axis by each function's parity ($\phi$ odd, the rest even) and through the equator evenly.

`field_equations` evaluates the paper's seven equations on the authors' grid by fourth order differences.
For $0.15 < \mu r < 30$, away from the axis, the largest residual is $1.5\times10^{-4}$ of the largest term of its equation, and $10^{-5}$ for the wave equation, which agrees with the error below $10^{-3}$ the paper's section 3.3 states.

From the data:

| | |
| --- | --- |
| mass, from $F_0 = c_t/r + \dots$ on the equator | $0.9330\,c^2/G\mu$ (the authors' page: $0.933$) |
| angular momentum, from $W = 2GJ/c^3r^3 + \dots$ | $0.7403\,c^3/G\mu^2$ ($0.739$) |
| horizon angular velocity, $cW(r_H)$ | $0.82\,\mu c$, constant on the horizon to $10^{-12}$ |
| horizon temperature, $e^{F_0 - F_1}/4\pi r_H$ | $0.02548\,\mu$ |
| horizon area | $3.647/\mu^2$ |
| Noether charge $Q$ and field energy $M_\Psi$ | $0.6258$ and $0.6987$ |
| horizon angular momentum $J - mQ$ and mass $M - M_\Psi$ | $0.1146$ and $0.2344$ ($0.114$ and $0.234$) |
| ergosurface on the equator | $r = 0.365/\mu$ |
| largest field, on the equator | $\sigma = 0.150$ at $r = 0.96/\mu$ |
| $e^{F_1 - F_0}$ on the horizon at the axis | $31.23$ |
| circumference radius of the horizon on the equator, $e^{F_2}r_H$ | $0.687/\mu$ |

Smarr's relation, $M = 2T_HS + 2\Omega_H(J - mQ) + M_\Psi$ with $S = A_H/4$, closes to $2\times10^{-6}$.

Two formulas of the paper of 2015 are misprinted, and the letter of 2014 has them right.
The temperature is $e^{F_0^{(0)} - F_1^{(0)}}/4\pi r_H$ with the horizon values $F^{(0)}$, the letter's form, where the paper's (3.12) prints $F^{(2)}$.
The field's energy (3.16) needs $e^{-2F_0}$ in $\mu^2 - 2e^{-2F_0}w(w - mW)/N$, where the paper prints $e^{-2F_2}$; with $e^{-2F_2}$ the field's energy comes out $-0.110$ and Smarr's relation fails by $0.81$.

Kerr's black hole of the same mass and angular momentum has $r_H = 0.98173/\mu$ and $c_t = -0.44217/\mu$, so $b = 0.44217/\mu$, `kerr_for`, a horizon area $33.39/\mu^2$, nine times the hairy hole's, and a horizon circumference radius $2GM/c^2 = 1.866/\mu$ on the equator.
Configuration IV lies outside the band between extremal Kerr and the existence line of the clouds, where the paper finds the hairy hole always the larger in area.

## Step 4: the drawings

- Spacetime diagrams, four views: the axis $\theta = 0$, where the dragging drops out and the rays are null geodesics, and the equator with $\varphi$ divided out, the shadows of the null geodesics with no angular momentum, the ergosurface and the radius of the largest field marked, each for configuration IV and for the Kerr member of the same mass and angular momentum. The rows declare the functions of each plane, `ksh_f0a` to `ksh_we`, from `Hair`; across either plane every function is even, so the rows hold the plane's metric and its first derivatives, and none reads the Kretschmann scalar.
- Conformal diagram, two views: the half axis outside the horizon of each, a diamond with the horizon on its left, through $p, q = \arctan((ct \mp r_*)/\ell)$. For the hair $r_*$ is `Hair.tortoise_fast`, $2r_He^{F_1 - F_0}\ln x$ plus a regular quadrature in $x$; for Kerr it is the Boyer-Lindquist tortoise coordinate of the axis in $R = r + b$.
- Embedding diagram, one view, turnable: the equatorial plane at one moment, $\rho = e^{F_2}r$ and $g_{rr} = e^{2F_1}/N$, through the bifurcation sphere into the mirror exterior, since the functions are even in $x$; Kerr's surface of the same mass and angular momentum is drawn beside it down to its own throat, set level with the hairy surface where both are $6/\mu$ in circumference radius. The solution is stationary, so there is no stack and no movie.

Only the outside of the horizon is computed, so no drawing goes inside it.

## Step 5: dimensions

`DIMENSIONS` declares $t$ a time, $r$, $r_H$ and $b$ lengths, $\theta$, $\varphi$ and the $F_i$ pure numbers, and $W$ an inverse length.
