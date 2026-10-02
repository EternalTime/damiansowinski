# Cremmer-Scherk spontaneous compactification

Minkowski space of four dimensions times a sphere of radius $a$, in six dimensions, held by a monopole's gauge field on the sphere and a cosmological constant.
The chart is written by `_tools/derivations/print_charts.py --metric cremmer_scherk`, which took under a second on 2 October 2026, and `verify_metrics.py --system cremmer_scherk/cartesian` checks it in under a second.

## Step 1. The papers and the source of the chart

Both papers were verified on Crossref and INSPIRE before anything was built: E. Cremmer and J. Scherk, Nucl. Phys. B 108, 409 (1976), doi:10.1016/0550-3213(76)90286-8, and Nucl. Phys. B 118, 61 (1977), doi:10.1016/0550-3213(77)90363-7.
Neither is open, and neither has a scan in the KEK library, so only their abstracts were read.
The line element and the relations between the radius, the couplings, and the cosmological constant are taken from Z. Horváth, L. Palla, E. Cremmer, and J. Scherk, Nucl. Phys. B 127, 57 (1977), read in full from the KEK scan of the preprint LPTENS 77/10 (`lib-extopc.kek.jp/preprints/PDF/1977/7706/7706154.pdf`).
That paper restates Cremmer and Scherk's solution in the Abelian gauge of Wu and Yang and says so: its results "correspond (with $G = SO(3)$ and $H = I_3$) to the particular case $n_i = (0, 0, 1)$" of the 1976 paper.

- `cartesian`: their (2), $ds^2 = -dx_0^2 + \sum_i dx_i^2 + R_0^2(d\theta^2 + \sin^2\theta\,d\varphi^2)$, with $R_0$ written $a$ and the flat coordinates written $t$, $x$, $y$, $z$.

No other chart of the six dimensions appears in the papers read.
The two regions of Wu and Yang are two gauges of the field over one chart of the sphere, so they are no second chart of the metric.
The 1977 paper's spheres of more dimensions are not published as charts: its text was not read, so their radii could not be sourced.

## Step 2. Their (4) and the cosmological constant

Their action is $-\frac{1}{16\pi G}\int\sqrt{-g}\,R - \int\sqrt{-g}\left[\frac{1}{4}\,2\,\mathrm{Tr}(G_{\mu\nu}G^{\mu\nu}) + \frac{1}{2}V_0\right]$, in units with $\hbar = c = 1$, and the field is $W_\varphi = (H/e)(\cos\theta \mp 1)$ on the two caps.
The field strength is $G_{\theta\varphi} = -(H/e)\sin\theta$, of the one magnitude $B = \sqrt{k}/(eR_0^2)$ in an orthonormal frame, with $k = 2\,\mathrm{Tr}H^2$.
Their (4) is $V_0 = e^2/(k(8\pi G)^2)$ and $R_0^2 = 8\pi Gk/e^2$.
The second is $8\pi GB^2 = 1/R_0^2$.
The term $\frac{1}{2}V_0$ is an energy density, so in $G_{\mu\nu} + \Lambda g_{\mu\nu} = 8\pi GT_{\mu\nu}$ it is $\Lambda = 8\pi GV_0/2 = e^2/(16\pi Gk) = 1/2R_0^2$.
The page states $\Lambda = 1/2a^2$ in that convention, which is also Plebański and Hacyan's $\Lambda = 1/2b^2$ for the flat plane times a sphere in four dimensions.

## Step 3. What the script checks

`cremmer_scherk_check` builds $F_{\theta\phi} = Ba^2\sin\theta$ and holds Maxwell's equations $\partial_\mu(\sqrt{-g}F^{\mu\nu}) = 0$, the stress $T^\mu{}_\nu = (B^2/2)\,\mathrm{diag}(-1, -1, -1, -1, 1, 1)$, and $G^\mu{}_\nu + \Lambda\delta^\mu{}_\nu = 8\pi GT^\mu{}_\nu$ at $8\pi GB^2 = 1/a^2$ and $\Lambda = 1/2a^2$.
The Einstein tensor is then $-1/a^2$ on the four flat dimensions and zero on the sphere.
An Abelian field stands for the monopole: the field of the paper lies along one generator $H$, and its stress is that of a Maxwell field of strength $B$.

## Step 4. The drawings

Everything is drawn in units of $a$.
The spacetime diagrams are the plane of $t$ and $x$, flat, and the plane of $t$ and $\phi$ on the equator of the sphere at one point of the flat dimensions, where $\phi = 0$ and $2\pi$ are one line and a ray returns after $2\pi a/c$.
The conformal diagram is the diamond of the plane of $t$ and $x$, as the plane $y = z = 0$ of `minkowski` is.
The embedding diagram has the cylinder of $x$ and $\phi$, Duff's hosepipe with a flat dimension for its length, and the sphere at one event.
The spacetime is static and homogeneous, so there is no movie and no stack.
The mass $\hbar\sqrt{J(J + 1)}/ac$ in the caption of the second spacetime diagram is the paper's $(J(J + 1) - q_m^2)^{1/2}/R_0$ at $q_m = 0$.
