# The quantum Oppenheimer-Snyder black hole

Jerzy Lewandowski, Yongge Ma, Jinsong Yang, and Cong Zhang's ball of dust in loop quantum cosmology, Physical Review Letters 130, 101501 (2023), arXiv:2210.02253, with the vacuum outside it,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2\,d\Omega^2,\qquad f = 1 - \frac{2m}{r} + \frac{\alpha m^2}{r^4},$$

with $m = GM/c^2$ and the area $\alpha = 16\sqrt{3}\,\pi\gamma^3\ell_p^2$.
Its five charts are written by `_tools/derivations/print_charts.py --metric quantum_oppenheimer_snyder`, and `verify_metrics.py --system quantum_oppenheimer_snyder/<chart>` checks each in seconds.
Equation numbers are those of the arXiv versions.

## Step 1. The parameters and the charts

Lewandowski, Ma, Yang, and Zhang write $GM$ and $\alpha G^2M^2$ with $c = 1$; with $m = GM/c^2$ no factor of $G$ or $c$ stands in a component, and $\alpha m^2/r^4$ is a pure number.
Each chart has a source.

- `static`, their (4): the static time $t$ and the areal radius $r$, determined by the junction for $r \ge r_b = (\alpha m/2)^{1/3}$, their (5).
- `painleve_gullstrand`: Jarod George Kelly, Robert Santacruz, and Edward Wilson-Ewing's effective metric in Painlevé-Gullstrand coordinates, Physical Review D 102, 106024 (2020), arXiv:2006.09302, Sec. IV B, valid for $x \ge x_{\min} = (\gamma^2\Delta R_S)^{1/3}$. With $R_S = 2m$ and $\gamma^2\Delta = \alpha/4$ their $F$ is our $f$ and $x_{\min}$ is $r_b$; both papers say the two metrics coincide.
- `eddington_finkelstein_ingoing` and `eddington_finkelstein_outgoing`: the charts on the tortoise coordinate, $dr_*/dr = 1/f$ with $r_* = 0$ at $r_b$. Lewandowski, Ma, Yang, and Zhang note that one advanced Eddington-Finkelstein coordinate covers both horizons the surface crosses on its way in.
- `interior_comoving`: their (1), the flat Friedmann metric of Ashtekar, Pawlowski, and Singh, with the scale factor left free in the chart and its law in the description.

## Step 2. What `quantum_os_check` holds

Outside, every chart carries the mixed Einstein tensor $G^t{}_t = G^r{}_r = -3\alpha m^2/r^6$ and $G^\theta{}_\theta = G^\phi{}_\phi = 6\alpha m^2/r^6$.
The first is their energy density $\rho^q = 3\alpha GM^2/8\pi r^6$; the trace gives the Ricci scalar $R = -6\alpha m^2/r^6$, and the Kretschmann scalar is printed in Kelly, Santacruz, and Wilson-Ewing's form, $K = (48m^2/r^6)(1 - 5\alpha m/r^3 + 39\alpha^2m^2/4r^6)$.
The static chart is Schwarzschild's published metric at $\alpha = 0$; $f = 1$ at $r_b$ and $f'(r_b)\,r_b = -6m/r_b$; $f$ and $f'$ vanish together at $r = 3m/2$ when $m^2 = 16\alpha/27$, which is their $M_{\min} = 4\sqrt{\alpha}/3\sqrt{3}\,G$, and the merged horizon is at $2\sqrt{\alpha/3}$; and in their $\beta$, $m^2 = 4\beta^4\alpha/(1 - \beta^2)^3$, $f$ vanishes at their $r_\pm$.
Each other exterior chart is the static one pulled back, along $c\,dt = c\,d\tau - N\,dr/f$ with $N^2 = 1 - f$ for Painlevé and Gullstrand's, and along $c\,dt = dv - dr/f$ and $du + dr/f$ for the others.

## Step 3. The interior and the junction

Their (2), $H^2 = (8\pi G/3)\rho(1 - \rho/\rho_c)$ with $\rho \propto a^{-3}$, has the solution Kelly, Santacruz, and Wilson-Ewing give for the radius of their star in Classical and Quantum Gravity 38, 04LT01 (2021), arXiv:2006.09325,

$$L(t) = \left[\gamma^2\Delta R_S\left(\frac{9t^2}{4\gamma^2\Delta} + 1\right)\right]^{1/3},$$

which in our letters is $R^3 = r_b^3 + 9mc^2\tau^2/2$, and so $a^3 = 1 + 9c^2\tau^2/\alpha$ with $a = 1$ at the bounce.
Then $8\pi G\rho/c^2 = 12/\alpha a^3$, and $\rho_c = 3c^2/2\pi G\alpha$ is their $\sqrt{3}/(32\pi^2\gamma^3G^2\hbar)$.
`quantum_os_check` holds the published $G^\tau{}_\tau$ with that $a$ to $-\kappa\rho(1 - \rho/\rho_c)$, no flux, and a vanishing Weyl tensor.
The surface $\chi_0^3 = \alpha m/2$ has $R = a\chi_0$, and $(dR/d(c\tau))^2 = 2m/R - \alpha m^2/R^4 = 1 - f(R)$: the radial geodesic of the exterior with $E = f\,dt/d\tau = 1$, which with the areal radius continuous is the junction of their appendix.

## Step 4. The surface in the drawings

`quantum_os.py` carries the surface into every exterior chart.
On the way in $dR/d\tau = -N$, the Painlevé-Gullstrand $T$ is $\tau$ and $dv/d\tau = 1/(1 + N)$; on the way out $dR/d\tau = +N$, $dv/d\tau = (1 + N)/f$ and $dT/d\tau = (1 + N^2)/f$, which run off to infinity as the surface reaches $r_-$ at $c\tau_- = 0.424\,m$ through the horizon of the white hole.
Near $\tau_-$ both rates have the pole $1/\kappa_-(\tau_- - \tau)$, since $N(r_-) = 1$, which is integrated in closed form.
The outgoing chart is the time reverse of the ingoing one and $R$ is even in $\tau$, so its surface is $R(u) = R_{\rm in}(-u)$.
The conformal time of the dust is $\eta = c\tau\,{}_2F_1(1/3, 1/2; 3/2; -9c^2\tau^2/\alpha)$.

The drawings take $m = 1$ and $\alpha = 5m^2/4$: $r_b = 0.855$, $r_- = 1.127$ and $r_+ = 1.777\,m$, $\kappa_-/\kappa_+ = 3.34$, and $1.16$ times the least mass.
At $\alpha = m^2$ the inner horizon would sit at exactly $r = m$, but $\kappa_-/\kappa_+ = 4.99$ squeezes the region inside $r_-$ into a corner of the conformal diagram, as Reissner-Nordström's does at $r_q = 0.4\,r_s$.

## Step 5. The conformal diagram of the collapse

The exterior keeps the tower's map.
Inside, the light rays are $\eta \mp \chi$ constant, and every one meets the surface, so a point of the ball is drawn at the $p$ the surface has where its outgoing ray arrives, $\eta - \chi + \chi_0$, and at the $q$ the surface has where its ingoing ray left, $\eta + \chi - \chi_0$.
Both families are then lines of constant $p$ or $q$, and `conformal.py` checks the map against the published interior metric like every other.
The ball is the lens between the surface and the centre, both running from $i^-$ of the lower exterior to $i^+$ of the upper one, as Lewandowski, Ma, Yang, and Zhang's Fig. 2(a) has it.
