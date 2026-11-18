# Gravitational waves in Friedmann universes

The code is `waves_in_friedmann.py` beside this file; the page is `MFS/assets/data/metrics/waves_in_friedmann.json`.

## Sources

- Bičák and Griffiths, Ann. Phys. 252, 180 (1996), arXiv gr-qc/9610041: the metric form (1.1), the five stiff backgrounds (2.4), (2.5), (3.15), (3.16), (4.7) with the scale factors (5.2) and the metric functions (5.4)-(5.8), the fluid potentials (5.11)-(5.13), the field equations (6.1)-(6.7), and the fronts of Secs. 7-11: $z = \eta$ (7.1), $\rho = \eta$ (8.1), $\mu = \eta$ (9.1), $\rho = \eta$ (10.1), $\zeta = \eta$ (11.x).
- Bičák and Griffiths, Phys. Rev. D 49, 900 (1994): the hyperboloid coordinates (2.6)-(2.13) of the open planar chart and the Friedmann normalisation (3.1)-(3.4).
- Alekseev and Griffiths, Phys. Rev. D 52, 4497 (1995): the profiles $H_n$ of (4.1)-(4.3), normalised to $H_n(-1) = 1$, and the closed form of $S$ for one term.

## Charts

Every line element is the source's multiplied by $-1$, for the signature $(-,+,+,+)$, with $\gamma$ kept as the area that sets the scale.
The radial coordinate of the cylindrical charts is written $r$ for the source's $\rho$, which the History uses for the energy density.
The open cylindrical chart writes $\tfrac12\sinh2r\tanh r = \sinh^2r$ and $\tfrac12\sinh2r\coth r = \cosh^2r$.

The plane wave chart is the region behind the front $z = \eta$ of a plane wave in the flat universe and ahead of the null plane $z = -\eta$, where Alekseev and Griffiths's $H_2(\zeta)$, written with $\arccos\zeta$ and $\sqrt{1 - \zeta^2}$, is real.
One term with $n = 2$ is taken: $V = A\eta^2H$, which solves $V_{\eta\eta} + V_\eta/\eta - V_{zz} = 0$, Bičák and Griffiths's (6.2) with $e^{-U} = 2\eta$.
For it (6.7) integrates to $S = -\tfrac18A^2\eta^2(4\eta^2H^2 + (\eta^2 - z^2)K^2)$ with $K = dH_2/d\zeta$, and $e^{-M} = 2\eta e^{-S}$; $H_2$ and $K$ vanish on the front, so $V$, $S$ and the metric are continuous across it, as the junction conditions of BG96 Sec. 6 require.
By BG96 (7.7) the lowest term $k = 2$ makes the Weyl tensor continuous on the front.
The tensors are printed with $H$, $K$, $V$ and $S$ held as names, their slopes in `RATES`: $H' = K$ and, from the profile's equation, $(1 - \zeta^2)K' + 3\zeta K - 4H = 0$.
The cylindrical and toroidal wave regions have no chart: their metric function $M$ is given in the sources only as an integral over the spectral amplitude.

Bičák and Griffiths's normalising constant $c_k$ of (6.11) is not used; Alekseev and Griffiths, Class. Quantum Grav. 14, 2869 (1997), footnote 1, write that it replaces incorrect expressions, and the profile here is normalised by the recursion and checked numerically.

## The fluid

Every chart is checked at points against $R_{ab} = \tfrac12\tilde\sigma_a\tilde\sigma_b$ with $\tilde\sigma = \sqrt3\log\eta$, $\sqrt3\log\tanh\eta$ and $\sqrt3\log\tan\eta$, the potentials of BG96 (5.11)-(5.13).
The factor $\tfrac12$ is the one the Friedmann normalisation of BG96 (5.1) and BG94 (3.1)-(3.4) fixes, $\rho R^6 = 3\gamma^2/8\pi$; BG96 (5.9) prints $16\pi\tilde\rho$ for what that normalisation makes $32\pi\tilde\rho$.
So $\rho = p = 3e^S/(64\pi\gamma\eta^3)$ in the flat universe, with $S = 0$ ahead of the front, $3/(8\pi\gamma\sinh^32\eta)$ in the open one and $3/(8\pi\gamma\sin^32\eta)$ in the closed one, and $R = -16\pi\rho$; the fluid is at rest in every chart.
`_tools/test_waves_in_friedmann.py` holds every published Einstein tensor to these densities.

## Drawings

Every plane of $\eta$ and one spatial coordinate carries a conformal factor times $-d\eta^2 + dx^2$, so the light rays are the lines $\eta \mp x$ constant in every chart and each wave front is the ray $x = \eta$.
The conformal diagrams are Minkowski's upper half plane, or its quarter $r \ge 0$, under $p, q = \arctan(\eta \mp x)$; the closed universe's square $\eta, \zeta \in (0, \pi/2)$ is drawn as it stands.
The embedding diagram is the closed universe at its greatest expansion, $\eta = \pi/4$, cut at $\sigma = 0$: $a^2(d\zeta^2 + \sin^2\zeta\,d\delta^2)$ with $a^2 = \gamma\sin2\eta = \gamma$, the hemisphere from the axis $\zeta = 0$ at its pole to the axis $\zeta = \pi/2$ at its rim, the fronts $\zeta$ constant its circles of latitude.
The curvature figure draws $K$ and $R$ across the wave region at $\eta = 1$ and $A = 2$ beside the flat universe's $K = 15/(16\gamma^2\eta^6)$ and $R = -3/(4\gamma\eta^3)$.
