# The black holes of string theory with three and four charges

The black holes whose entropy string theory first counted: in five dimensions with three charges and in four with four, each charge entering through a harmonic function of its own.
`print_charts.string_bh_charges` writes the five charts, `string_bh_charges_check` holds each to its field equations and to the published charts it reduces to, and `_tools/test_string_bh_three_four_charges.py` holds the published files and the drawings to their numbers.

## Step 1. The sources, each verified before anything was built

Every paper of the candidate entry was found at Crossref and at INSPIRE, and each line element below was read off the arXiv copy.

- Andrew Strominger and Cumrun Vafa, Phys. Lett. B 379, 99 (1996), doi:10.1016/0370-2693(96)00345-0, arXiv:hep-th/9601029. Their (2.8) is the extreme Reissner-Nordström hole of five dimensions in the areal radius.
- Curtis G. Callan and Juan M. Maldacena, Nucl. Phys. B 472, 591 (1996), doi:10.1016/0550-3213(96)00225-8, arXiv:hep-th/9602043. Their (2.1) is the Reissner-Nordström hole of five dimensions, $\lambda = (1 - r_+^2/r^2)(1 - r_-^2/r^2)$, and (2.11) the extreme hole of three charges.
- Gary T. Horowitz, Juan M. Maldacena and Andrew Strominger, Phys. Lett. B 383, 151 (1996), doi:10.1016/0370-2693(96)00738-1, arXiv:hep-th/9603109. Their (2.7) and (2.8) are the chart of three charges.
- Mirjam Cvetič and Donam Youm, Phys. Rev. D 53, R584 (1996), doi:10.1103/PhysRevD.53.R584, arXiv:hep-th/9507090. Their (9) is the extreme chart of four charges.
- A. A. Tseytlin, Mod. Phys. Lett. A 11, 689 (1996), doi:10.1142/S0217732396000709, arXiv:hep-th/9601177. His (31) is the extreme chart of three charges.
- A. A. Tseytlin, Nucl. Phys. B 475, 149 (1996), doi:10.1016/0550-3213(96)00328-8, arXiv:hep-th/9604035. Crossref misprints the name as Tseytfin; INSPIRE and the arXiv copy have Tseytlin.
- Juan M. Maldacena and Andrew Strominger, Phys. Rev. Lett. 77, 428 (1996), doi:10.1103/PhysRevLett.77.428, arXiv:hep-th/9603060.
- Clifford V. Johnson, Ramzi R. Khuri and Robert C. Myers, Phys. Lett. B 378, 78 (1996), doi:10.1016/0370-2693(96)00383-8, arXiv:hep-th/9603061.
- Gary T. Horowitz, David A. Lowe and Juan M. Maldacena, Phys. Rev. Lett. 77, 430 (1996), doi:10.1103/PhysRevLett.77.430, arXiv:hep-th/9603195. Their (2) is the chart of four charges, after Cvetič and Youm, Nucl. Phys. B 472, 249 (1996), arXiv:hep-th/9512127, whose list of special cases gives the dilaton couplings $a = 0$, $1/\sqrt3$, $1$ and $\sqrt3$ for four, three, two and one equal charges.

## Step 2. The charts

Horowitz, Maldacena and Strominger write boosts $\alpha$, $\gamma$ and $\sigma$; the charts write the lengths $r_i = r_0\sinh\alpha_i$ in five dimensions and $r_i = r_0\sinh^2\alpha_i$ in four, so that $H_i = 1 + r_i^2/r^2$ and $H_i = 1 + r_i/r$ and the extreme hole is $r_0 = 0$ with the $r_i$ kept.

- `five_charges`: $-(H_1H_2H_3)^{-2/3}f\,c^2dt^2 + (H_1H_2H_3)^{1/3}(dr^2/f + r^2d\Omega_3^2)$, $f = 1 - r_0^2/r^2$.
- `five_extreme`: the same at $f = 1$.
- `five_areal`: the three charges equal, $r_i = r_q$, in $\rho^2 = r^2 + r_q^2$, where $r_+^2 = r_0^2 + r_q^2$ and $r_-^2 = r_q^2$.
- `four_charges`: $-f\,c^2dt^2/\sqrt{H_1H_2H_3H_4} + \sqrt{H_1H_2H_3H_4}\,(dr^2/f + r^2d\Omega^2)$, $f = 1 - r_0/r$.
- `four_extreme`: the same at $f = 1$.

A chart that names $f$ and the $H_i$ holds them as functions, `HELD` in `verify_metrics.py`, with the slopes `RATES` declares: $\partial_rH_i = -2(H_i - 1)/r$ and $\partial_rf = 2(1 - f)/r$ in five dimensions, and half of each in four.
Every value is then a rational function of $r$, $f$ and the $H_i$ with a root of the $H_i$ in front.
Written out in $r$ and the $r_i$, the chart of three charges took 151 seconds and printed a Riemann component as a polynomial of fifty terms; held, it takes seven seconds.

## Step 3. The field equations

Each chart solves $R - \tfrac12(\partial\phi)^2 - \tfrac14\sum_iX_i^{-2}F_i^2$ with $n = D - 2$ gauge fields and scalars $X_i = H_i^{-1}(H_1\cdots H_n)^{1/n}$, whose product is 1: $A^i = q_i\,dt/(r^{D-3}H_i)$ with $q_i^2 = r_i^2(r_i^2 + r_0^2)$ in five dimensions and $r_i(r_i + r_0)$ in four.
`string_bh_charges_check` checks Einstein's equations slot by slot, Maxwell's equation of each field and the equation of each scalar, at three random points to thirty digits, each held name standing as a number of its own.
With the charges equal every $X_i$ is 1 and the equations are the Einstein-Maxwell equations, which the areal chart is checked against directly.

## Step 4. The charts it reduces to

The same check compares each chart with a published one: Tangherlini's and Schwarzschild's with no charge; `rn_metric` in $r + r_q$ with the four charges equal; `dilaton_black_hole` in $r + r_q$ with two charges equal and two absent; the Einstein metric of `kaluza_klein_black_hole` with one charge; and the single hole of `majumdar_papapetrou` for the extreme chart with the four charges equal.
The area of the horizon is $2\pi^2\prod\sqrt{r_0^2 + r_i^2}$ in five dimensions and $4\pi\prod\sqrt{r_0 + r_i}$ in four, Horowitz, Maldacena and Strominger's (2.9) and Horowitz, Lowe and Maldacena's (5), and $2\pi^2r_1r_2r_3$ and $4\pi\sqrt{r_1r_2r_3r_4}$ at $r_0 = 0$.

## Step 5. The drawings

`slices.SBC` names the parameters every drawing uses: $r_0 = 1$ with $r_i = 1/2, 1, 3/2$ and, in four dimensions, $2$; the extreme holes in units of $r_2$ with $(1/2, 1, 2)$ and $(1/2, 1, 3/2, 2)$; and the areal chart at $r_0 = 3r_q/4$, so that its horizons are $5r_q/4$ and $r_q$.
The tortoise coordinate of a chart with unequal charges is an elliptic integral, so `slices.sbc_rstar` writes it as its poles plus a smooth rest: with $S$ the root of $\prod(r^2 + r_i^2)$ or of $\prod(r + r_i)$ and $L$ the polynomial that carries the poles, the rest is $(S^2 - L^2)/(D(S + L))$, where $S^2 - L^2$ is a polynomial the poles divide exactly, and it is integrated by quadrature.
The spacetime diagrams hold every ray to $ct \mp r_*$, the conformal diagrams draw Kruskal's block of each chart off extremality with the inner horizon $r = 0$ for its upper edges, the whole Reissner-Nordström tower for the areal chart, and the exterior diamond of each extreme hole, and the embedding diagrams draw the moment $t = 0$ of each chart.
