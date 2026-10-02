# Bowers and Liang's anisotropic star: the chart, the source and what the diagrams draw

Bowers and Liang's star is a static sphere of uniform density whose pressure across the radius differs from its pressure along the radius.
This note records where the chart of `bowers_liang.json` comes from, how the checker reads its two held names, and what the diagrams take from it.
Units here have $G = c = 1$ unless a line keeps them; the file keeps both.

## Step 1. The equations of the paper

Bowers and Liang (1974) write a static sphere as $ds^2 = e^\nu dt^2 - e^\lambda dr^2 - r^2d\theta^2 - r^2\sin^2\theta\,d\varphi^2$, their (2.2), with the stresses $T^\mu{}_\nu = \mathrm{diag}(\rho, -p_r, -p_\perp, -p_\perp)$, their (2.3).
Their (2.7) is the equation of equilibrium, (2.8) the slope of $\nu$, and (2.10) the radial component,

$$\frac{dp_r}{dr} = -(\rho + p_r)\frac{\nu'}{2} + \frac{2}{r}(p_\perp - p_r),
\qquad \frac{\nu'}{2} = \frac{m + 4\pi r^3p_r}{r(r - 2m)},
\qquad e^{-\lambda} = 1 - \frac{2m}{r},$$

with $m(r) = \int_0^r 4\pi r^2\rho\,dr$, their (2.9).
For the incompressible star, $\rho = \rho_0$, they take their (3.1) with $n = 2$ and (3.2),

$$p_\perp - p_r = C\,\frac{r^2(\rho_0 + p_r)(\rho_0 + 3p_r)}{1 - 2m/r},$$

with $C$ a constant, which turns (2.7) into their (3.3) and integrates to their (3.4),

$$p_r = \rho_0\,\frac{(1 - 2m/r)^Q - (1 - 2M/R)^Q}{3(1 - 2M/R)^Q - (1 - 2m/r)^Q}, \qquad Q = \frac{1}{2} - \frac{3C}{4\pi}.$$

The paper was read in full from the scan the Astrophysics Data System holds.

## Step 2. The chart

The paper stops at the pressure and never writes $e^\nu$.
With $m = Mr^3/R^3$ the file writes $Z = e^{-\lambda} = 1 - r_sr^2/R^3$, with $r_s = 2M$, and $Z_R = 1 - r_s/R$ for its value at the surface.
Then $8\pi\rho_0 = 3r_s/R^3$, and (2.8) with (3.4) is

$$\frac{\nu'}{2} = \frac{r_sr}{R^3}\,\frac{Z^{Q-1}}{3Z_R^Q - Z^Q} = \frac{1}{2Q}\,\frac{d}{dr}\ln\left(3Z_R^Q - Z^Q\right),$$

so $e^\nu$ is a constant times $(3Z_R^Q - Z^Q)^{1/Q}$, and the constant is the one that makes $e^\nu = Z_R$ at the surface:

$$e^\nu = N^{1/Q}, \qquad N = \frac{1}{2}\left(3Z_R^Q - Z^Q\right).$$

At $Q = 1/2$ this is Schwarzschild's $\tfrac{1}{4}(3\sqrt{Z_R} - \sqrt{Z})^2$.
Dev and Gleiser (2003) write the same function in their (101), with the factor $\tfrac{1}{4}$ in front of $(3y_1^{2Q} - y^{2Q})^{1/Q}$ where the surface asks for $2^{-1/Q}$; the two agree at $Q = 1/2$ alone.
The file's third name is

$$P = \frac{Z^Q}{N} = 1 + \frac{3p_r}{\rho_0},$$

by (3.4), which is 1 at the surface.
`print_charts.bowers_liang` builds the chart, and `bowers_liang_check` holds it to the uniform density, to (3.4), to (3.1) with (3.2), and to Schwarzschild's published exterior at $r = R$.

## Step 3. Two names held with their slopes

$N$ and $P$ hold a power $Q$ of $Z$, and each has a slope that is algebraic in the names themselves:

$$\partial_rN = \frac{Qr_sr\,NP}{R^3Z}, \qquad \partial_rP = -\frac{Qr_sr\,P(P + 2)}{R^3Z}.$$

`HELD` and `RATES` in `verify_metrics.py` declare both, the reader checks each slope against the definition the file publishes, and from there every tensor is a rational function of $r$, $R$, $r_s$, $Q$ and $P$, with one factor $N^{1/Q}$ for each lowered time index.
$P$ is defined in the file by $Z$ alone, $P = 2Z^Q/(3Z_R^Q - Z^Q)$, since the drawings write a held name's definition in once, and a definition that named $N$ would leave $N$ standing.
The printer writes $r_sr^2 = R^3(1 - Z)$ inside every sum and groups the sum by powers of $Z$ and of $1 - Z$.
Each term the anisotropy adds then carries $(1 - Z)(2Q - 1)(P + 2)$: it vanishes at $Q = 1/2$ and at the centre, and the Weyl tensor is that factor alone, $C^t{}_{rtr} = r_s^2r^2P(2Q - 1)(P + 2)/(12R^6Z^2)$.

A power $Q$ of $1 - r_sr^2/R^3$ put through the checker's canonical form is spelled $(-1)^Q(r_sr^2 - R^3)^Q/R^{3Q}$, which is sound in every comparison and is not the number the power has.
So the drawings and `_tools/test_bowers_liang.py` read the published components as they stand and write the held definitions in without the canonical form.

## Step 4. The limits of the family

The central pressure is infinite where $N$ vanishes at the centre, $3Z_R^Q = 1$, which is $r_s/R = 1 - 3^{-1/Q}$, their (3.6): $8/9$ at $Q = 1/2$, $80/81$ at $Q = 1/4$, and 1 as $Q \to 0$.
The surface redshift there is $3^{1/2Q} - 1$, their (3.12), and at one density the mass goes as $(r_s/R)^{3/2}$, so the heaviest star holds $(9/8)^{3/2} = 1.19$ of Schwarzschild's heaviest, their (3.10).

As $Q \to 0$, $N^{1/Q} \to Z_R^{3/2}Z^{-1/2}$ and $p_r \to 0$: the star is Florides's (1974) ball of uniform density with no stress along the radius, the uniform chart of `einstein_cluster.json`.
The check takes that limit of the chart and compares it with the published $g_{tt}$ there.

The paper writes its critical compactness by $\xi$, $(2M/R)_{\rm crit} = 1 - (1/3)^{2/(1 - \xi)}$, and prints $\xi \equiv 3C/2$.
With $Q = 1/2 - 3C/4\pi$, which its (3.3) confirms, $1/Q = 2/(1 - \xi)$ asks for $\xi = 3C/2\pi$.
Yagi and Yunes (2015) agree: their $\lambda_{\rm BL} = -3C$ reaches the black hole limit at $-2\pi$, which is $C = 2\pi/3$ and $Q = 0$.
The paper's "we need only $C = 0.063$" for the redshift 2.358 uses the printed $\xi = 3C/2$; in the exponent it is $Q = \ln 3/(2\ln 3.358) = 0.453$, which is what the History quotes.
The file takes $Q$ for its parameter, which the slip does not touch.

## Step 5. What the diagrams draw

Every diagram draws the star with $Q = 1/4$ at Schwarzschild's limit, $R = 9r_s/8$, where the star of equal pressures has an infinite central pressure.
There $Z_R^Q = 1/\sqrt{3}$, so the central pressure is $\rho_0/\sqrt{3}$, and at the surface $p_\perp = \rho_0$ exactly, since $p_\perp(R)/\rho_0 = (1 - 2Q)(r_s/R)/(4Z_R)$.
The stresses add up to more than the density, as Andréasson's (2008) bound asks of a star at $r_s/R = 8/9$: with $p_r + 2p_\perp \le \Omega\rho$ a static sphere has $r_s/R \le 1 - 1/(1 + 2\Omega)^2$.
$-g_{tt}$ runs from $0.0179$ at the centre to $1/9$ at the surface, light moves at $dr/d(ct) = \sqrt{-g_{tt}Z}$ between $0.134$ and $0.111$, and it crosses from the centre to the surface in $8.28\,r_s/c$ of $t$, which `conformal.bowers_liang` checks by quadrature.
The spacetime diagrams are drawn twice as tall as wide, the limit, to show the lean of those cones, and the conformal diagram scales $ct \mp r_*$ by $9r_s$ so that the star, $8.28\,r_s$ deep in $r_*$, takes a third of Minkowski's triangle.

On a slice of constant $t$ the equator has $g_{rr} = 1/Z$ for every $Q$, the metric of a sphere of radius $\sqrt{R^3/r_s}$, so the embedded surface is the cap Schwarzschild's star has, which `embedding.bowers_liang` checks against $z = a - \sqrt{a^2 - r^2}$.
The cap ends $\arcsin\sqrt{r_s/R} = 70.5°$ from its pole at $R = 9r_s/8$.
