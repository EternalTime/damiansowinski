# Bekenstein's universes with no singularity

Friedmann universes of dust whose particles carry the mass $\mu + f\psi$, a conformally coupled massless scalar field $\psi$ and radiation, in three charts, one for each sign of the spatial curvature.
`print_charts.py` writes every chart's mathematics and `bekenstein_universes_check` holds each to the field equations below at points.
`_tools/test_bekenstein_universes.py` holds the conditions and limits the page states.

## Sources

Bekenstein, Phys. Rev. D 11, 2072 (1975), is the founding paper.
Every equation here is derived from the action and checked against the restatement of Novello and Perez Bergliaffa, Phys. Rep. 463, 127 (2008), section 3.1, equations (64)-(67), which it matches in every coefficient.

## Step 1: the field equations

In units $G = c = 1$ and signature $(-,+,+,+)$ the action is
$$S = \int\sqrt{-g}\left(\frac{R}{16\pi} - \frac{1}{2}\nabla_\alpha\psi\nabla^\alpha\psi - \frac{R\psi^2}{12}\right)d^4x - \sum\int(\mu + f\psi)\,d\tau,$$
plus radiation.
The field's stress tensor is the conformal one,
$$T^{\psi}_{\mu\nu} = \nabla_\mu\psi\nabla_\nu\psi - \frac{1}{2}g_{\mu\nu}(\nabla\psi)^2 + \frac{1}{6}\left(g_{\mu\nu}\Box - \nabla_\mu\nabla_\nu + G_{\mu\nu}\right)\psi^2,$$
which is traceless on the field's own equation, $\Box\psi - R\psi/6 = fn$.
The dust is at rest with number density $n = N/a^3$ and energy density $n(\mu + f\psi)$, and the radiation has density $B/a^4$.

In the metric $ds^2 = a^2(\eta)\left(-d\eta^2 + d\chi^2 + S_k^2(\chi)d\Omega^2\right)$, with $F = a\psi$ and $\ell = 4\pi N\mu/3$:
$$F'' + kF = -fN,\qquad a'' + ka = \ell,\qquad a'^2 + ka^2 = \frac{4\pi}{3}\left(F'^2 + kF^2 + 2N\mu a + 2NfF + 2B\right).$$
The second is the trace, $R = 6(a'' + ka)/a^3 = 8\pi N\mu/a^3$: only the bare mass enters it.

## Step 2: the solutions

The time translation puts the bounce at $\eta = 0$ in each chart.

- Closed, $k = 1$: $a = \ell - A\cos\eta$ and $F = -fN + P\cos\eta + Q\sin\eta$, with $A^2 = \ell^2 + (4\pi/3)(2B + P^2 + Q^2 - f^2N^2)$. The radius stays positive for all $\eta$ when $A < \ell$, that is when $f^2N^2 > 2B + P^2 + Q^2$. Novello and Perez Bergliaffa write the condition $N^2f^2 > 2B$, the case $P = Q = 0$.
- Flat, $k = 0$: $a = a_b + \ell\eta^2/2$ and $F = -fN\eta^2/2 + d\eta + e$, with $2N\mu a_b = -(2B + d^2 + 2Nfe)$. The bounce radius $a_b$ is positive when $-2Nfe > 2B + d^2$.
- Open, $k = -1$: $a = A\cosh\eta - \ell$ and $F = fN + D_+e^\eta + D_-e^{-\eta}$, with $A^2 = \ell^2 - (4\pi/3)(2B + f^2N^2 - 4D_+D_-)$. The radius stays positive when $A > \ell$, that is when $4D_+D_- > 2B + f^2N^2$.

The total energy density, $3(a'^2 + ka^2)/8\pi a^4$, is at the bounce $3k/8\pi a^2$: positive closed, zero flat, negative open.

Limits: $A = 0$ is a static universe of radius $\ell$, the Einstein static geometry with no cosmological constant; $f = B = P = Q = 0$ gives $A = \ell$ and Friedmann's closed dust, $\ell(1 - \cos\eta)$; $a_b = 0$ gives flat dust, $\ell\eta^2/2$; and $A = \ell$ in the open chart gives open dust, $\ell(\cosh\eta - 1)$.

## Step 3: the charts

Each chart prints its components in the scale factor $a$, which `bekenstein_scale_factor` substitutes for the even part of the time dependence, leaving only the odd factor $\sin\eta$, $\eta$ or $\sinh\eta$, the rate of change of $a$.
The Ricci scalar is $6\ell/a^3$ in every chart and the Weyl tensor vanishes.
`bekenstein_universes_check` picks the constants $(f, P, Q, B)$ inside each chart's range with no singularity, fixes the amplitude from the Friedmann constraint, and compares $G_{\mu\nu}$ with $8\pi T_{\mu\nu}$, the field equation and $R = 8\pi N\mu/a^3$ at three points in thirty digits, with no symbolic determinant or simplification.

## Step 4: the drawings

The drawings take $\ell$ as their unit of length.
Closed: $A = 1/2$, so the radius runs between $1/2$ and $3/2$.
Flat: $a_b = 1/2$.
Open: $A = 3/2$, so the bounce radius is $1/2$.
