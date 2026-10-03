# Kerr's black hole in Bertotti and Robinson's field

Podolský and Ovcharenko's metric, the line element of `podolsky2025` (arXiv:2507.05199, read in its TeX source, its label `Kerr-BR`), is

$$ds^2 = \frac{1}{\Omega^2}\left[-\frac{Q}{\rho^2}\left(dt - a\sin^2\theta\,d\varphi\right)^2 + \frac{\rho^2}{Q}dr^2 + \frac{\rho^2}{P}d\theta^2 + \frac{P}{\rho^2}\sin^2\theta\left(a\,dt - \left(r^2 + a^2\right)d\varphi\right)^2\right],$$

with the metric functions that follow it,

$$\rho^2 = r^2 + a^2\cos^2\theta,\quad P = 1 + B^2\left(m^2\frac{I_2}{I_1^2} - a^2\right)\cos^2\theta,\quad Q = \left(1 + B^2r^2\right)\Delta,\quad \Omega^2 = 1 + B^2r^2 - B^2\Delta\cos^2\theta,$$

$$\Delta = \left(1 - B^2m^2\frac{I_2}{I_1^2}\right)r^2 - 2m\frac{I_2}{I_1}r + a^2,\qquad I_1 = 1 - \tfrac{1}{2}B^2a^2,\quad I_2 = 1 - B^2a^2.$$

Both charts are written by `_tools/derivations/print_charts.py --metric kerr_bertotti_robinson`, one chart to a run with `--system`; on 2 October 2026 the static chart took 7 seconds and the spinning chart 40, and `verify_metrics.py --system kerr_bertotti_robinson/static --system kerr_bertotti_robinson/boyer_lindquist` checks both in 15.

## Step 1. The period of the azimuth

Their angle $\varphi$ runs over $[0, 2\pi C)$ with the conicity of their section on the regularity of the axes, $C = 1/P(0) = 1/P(\pi)$, which makes both halves of the axis regular at once.
Both charts write $\varphi = C\phi$ with $\phi$ of period $2\pi$, as `kerr_melvin` writes Ernst and Wild's angle, so that every drawing reads a circumference from the published $g_{\phi\phi}$ as it does for Kerr.
$C$ is a constant of each chart and no defined name: written out, it stood as a quotient of polynomials in $a$, $m$ and $B$ in every component that holds $d\phi$.
Its value is stated in its description, and `kerr_bertotti_robinson_conicity` holds $CP$ to 1 at both poles.

## Step 2. The spinning chart

Held as Podolský and Ovcharenko write it, with $\rho^2$, $P$, $Q$ and $\Omega$ as functions, the Riemann tensor took a minute and ran to three quarters of a million operations.
Their section on the electromagnetic field gives the observers of no angular momentum, with the squared lapse $N^2 = PQ\rho^2/(R\,\Omega^2)$ and the rate $\omega = a(P(r^2 + a^2) - Q)/R$, $R = P(r^2 + a^2)^2 - Qa^2\sin^2\theta$.
The chart writes the line element in those four functions of a stationary field with an axis,

$$ds^2 = -N\,c^2dt^2 + F\left(\frac{dr^2}{Q} + \frac{d\theta^2}{P}\right) + W\left(C\,d\phi - \omega\,c\,dt\right)^2,\qquad N = \frac{PQ\Sigma}{A\Omega^2},\quad F = \frac{\Sigma}{\Omega^2},\quad W = \frac{A\sin^2\theta}{\Sigma\Omega^2},$$

with their $\rho^2$ written $\Sigma$ and their $R$ written $A$, so that neither collides with the radius or the curvature, and their $N^2$ written $N$, the squared lapse, as `kerr_melvin` writes it.
$Q$, $P$, $\omega$, $N$, $F$ and $W$ are held as functions while the tensors are built, `HELD` in `verify_metrics.py`, and the Riemann tensor then takes four seconds.
`kerr_bertotti_robinson_check` writes the names out and holds the chart to their line element slot by slot.
The Ricci scalar is published in the names, a sum of their second derivatives that no relation among six free functions makes zero; for the definitions it vanishes, the trace of Maxwell's stress, which the check holds at two points and `_tools/test_kerr_bertotti_robinson.py` at two more, in 30 digits.
The Kretschmann scalar is published in the names as well.
The checker compares each published value with the computed one as functions of the names, `agree_held`, which `kerr_melvin` brought in.

## Step 3. The static chart

At $a = 0$ their line element of the Schwarzschild-Bertotti-Robinson case is

$$ds^2 = \frac{1}{\Omega^2}\left(-f\,c^2dt^2 + \frac{dr^2}{f} + r^2\left(\frac{d\theta^2}{P} + C^2P\sin^2\theta\,d\phi^2\right)\right),\quad P = 1 + B^2m^2\cos^2\theta,\quad f = \left(1 + B^2r^2\right)\left(1 - B^2m^2 - \frac{2m}{r}\right),$$

with their calligraphic $\mathcal{Q}$ written $f$, since it is a pure number where the spinning chart's $Q$ is an area.
They print $\Omega^2 = 1 + B^2\left[r^2\sin^2\theta + \left(2mr + B^2m^2r^2\right)\cos^2\theta\right]$, which is the general $\Omega^2$ at $a = 0$; the chart writes it as $1 + B^2r^2 - B^2\left(\left(1 - B^2m^2\right)r^2 - 2mr\right)\cos^2\theta$.
$\Omega$, $f$ and $P$ are held, the Ricci scalar is stated as $0$, which the printer checks with the names written out, and the Kretschmann scalar is stated as

$$K = \frac{48m^2\left(1 + B^2mr\cos^2\theta\right)^2\Omega^4}{r^6} + \frac{8B^4\left(fP\sin^2\theta + \left(D_0 + D_1r\right)^2\cos^2\theta\right)^2}{\Omega^4},$$

$D_0 = 1 - B^2m^2(1 - 2\cos^2\theta)$ and $D_1 = B^2m(2 - (1 - B^2m^2)\cos^2\theta)$ being the coefficients of their Maxwell scalars at $a = 0$.
The first term is $48\Psi_2^2$ with their Weyl scalar at $a = 0$, $\Psi_2 = -(m/r^3)(1 + B^2mr\cos^2\theta)\Omega^2$, and the second is $2R_{ab}R^{ab} = 128|\Phi_0\Phi_2 - \Phi_1^2|^2$ with their Maxwell scalars; the constant 128 is fixed by the value $8B^4$ of Bertotti and Robinson's universe at $m = 0$, and the printer checks the whole against sympy.
The check holds the chart at $B = 0$ to the published spherical chart of `schwarzschild` with $r_s = 2m$, its horizon to $r_h = 2m/(1 - B^2m^2)$ with the surface gravity $(1 + B^2m^2)^2/(4m)$, both from their sections on horizons and thermodynamics, and its Weyl tensor to vanishing at $m = 0$.

## Step 4. The field

Their potential, at the duality angle $\gamma = 0$, is the complex one-form

$$\mathbf{A} = \frac{1}{2B}\left[\Omega_{,r}\,\frac{a\,dt - \left(r^2 + a^2\right)d\varphi}{r + ia\cos\theta} + \frac{i\,\Omega_{,\theta}}{\sin\theta}\,\frac{dt - a\sin^2\theta\,d\varphi}{r + ia\cos\theta} + \left(\Omega - 1\right)d\varphi\right],$$

and the field is $d$ of twice its real part.
`kerr_bertotti_robinson_einstein_maxwell` holds $G_{\mu\nu} = 2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}F^2\right)$ in geometric units at two points of each chart, to 30 digits relative to the largest component of $G_{\mu\nu}$; the real part is taken once the numbers are in, since `expand_complex` on the symbolic potential did not finish in five minutes.
At $a = 0$ the term in $i$ is purely imaginary, and the static chart's potential is $(C/B)(\Omega - 1 - r\,\partial_r\Omega)\,d\phi$.

## Step 5. The charge

Ovcharenko and Podolský showed in arXiv:2608.21672 (`ovcharenko2026b`, their section on the special charge) that the metric of `podolsky2025` is the member $e = e_s = \tilde m aB/\sqrt{1 - a^2B^2}$ of their Kerr-Newman-Bertotti-Robinson family, whose charge is $q_e = Ce$, with $\tilde m = m I_1/I_2$ their mass parameter.
In the $m$ of these charts that is $q_e = C\,maB\sqrt{I_2}/I_1$, which the spin parameter's description states in Gaussian units.
`_tools/test_kerr_bertotti_robinson.py` integrates the flux of the dual of the published field through spheres at $r = 2.5\,m$ and $7\,m$ and holds it to that value to twelve places, and the flux of the field itself to zero.
At the drawings' values it is $0.1405\,m$, against $2amB = 0.4\,m$ for Ernst and Wild's hole of the same $a$, $m$ and $B$.

## Step 6. The values the drawings are drawn at

The spinning hole is drawn at $m = 1$, $a = 4m/5$ and $B = 1/(4m)$, the values `kerr_melvin` draws Ernst and Wild's hole at, so that $I_1 = 49/50$, $I_2 = 24/25$ and $C = 60025/61374$.
Their formula for the roots of $\Delta$ puts the horizons at $r_+ = 1.6845\,m$ and $r_- = 0.4053\,m$, and on the equator $g_{tt}$ vanishes where $Q = a^2P$, their condition for the edge of the ergoregion, at $r = 2.0210\,m$.
The hole with no spin is drawn at $B = 1/(4m)$, where $C = 16/17$, $r_h = 32m/15$, and the innermost stable circular orbit is at $3r_h = 32m/5$, their result for the orbits of the static hole; the test finds it from the published equator as the least angular momentum of a circular orbit.

The embedded equator of the spinning hole has the circumference radius $C\sqrt{W}$, which on the equator is $C\sqrt{A}/(r\sqrt{1 + B^2r^2})$.
It is $1.8608\,m$ at $r_+$ and grows without a turn toward $C\sqrt{1 - a^2B^2\left(1 - B^2m^2I_2/I_1^2\right)}/B = 3.8380\,m$ at infinity, while $g_{rr} = F/Q$ falls as $r^{-4}$, so the radius $r$ runs to infinity within $7.5602\,m$ of the throat.
Their footnote on the Bertotti-Robinson case says that $r \to \infty$ is not the conformal infinity, since $\Psi_2$ does not vanish there, and that the conformal infinity and the causal structure were left for later work; no conformal diagram is drawn, and `conformal.py` lists the spacetime in `NOT_DRAWN`.
