# The Sultana-Dyer black hole

Sultana and Dyer's black hole in the Einstein-de Sitter universe is Schwarzschild's metric multiplied by the square of the scale factor of a flat universe of dust, $a = \eta^2/\eta_0^2$ in its conformal time,

$$ds^2 = \frac{\eta^4}{\eta_0^4}\left[-d\eta^2 + dr^2 + r^2d\Omega^2 + \frac{r_s}{r}\left(d\eta + dr\right)^2\right].$$

Its three charts are written by `_tools/derivations/print_charts.py --metric sultana_dyer`, and `verify_metrics.py --system sultana_dyer/<chart>` checks each in seconds.

## Step 1: the parameters and the units

Sultana and Dyer, Gen. Rel. Grav. 37, 1347, and Saida, Harada and Maeda, Class. Quantum Grav. 24, 4711, write the metric with $G = c = 1$, a mass $m$ or $M$ and $a = \eta^2$ or $(\eta/\eta_*)^2$.
The charts quote $r_s = 2GM/c^2$ for $2M$ and keep Saida, Harada and Maeda's constant as $\eta_0$, the conformal time at which $a = 1$.
The conformal time is a length, as the conformal chart of `frw` has it, with $a\,d\eta = c\,dt$ for the cosmic time, so no $c$ stands in the Kerr-Schild chart.
Faraoni, Phys. Rev. D 80, 044013, notes that Sultana and Dyer's own paper writes $t$ for the conformal time and $\eta$ for the cosmic time; the charts follow Saida, Harada and Maeda and Faraoni, with $\eta$ conformal.

## Step 2: the sources of the charts

- `kerr_schild`: Sultana and Dyer's own form, as Saida, Harada and Maeda's section 3.1 and Majhi's section 2, JCAP 05 (2014) 014, print it: flat space plus $(r_s/r)(d\eta + dr)^2$, all times $a^2$. Every slice of constant $\eta$ is spacelike, $g^{\eta\eta} = -(1 + r_s/r)/a^2$, and crosses the horizon.
- `schwarzschild_time`: $\eta = ct + r_s\ln(r/r_s - 1)$, the transformation Saida, Harada and Maeda and Majhi both write, which makes the metric $a^2$ times Schwarzschild's in its own coordinates. Carrera and Giulini, Phys. Rev. D 81, 043521, write the conformal factor in these coordinates as $(T + 2M_0\ln(R/2M_0 - 1))^2$. The chart names $\eta$ and `SultanaDyerForms` prints every value in it, the time entering through $\eta$ alone.
- `eddington_finkelstein_ingoing`: Majhi's advanced time $v = t + r_*$, with $r$ kept as the second coordinate as the ingoing chart of `schwarzschild` keeps it, so that $\eta = v - r$ and the conformal factor is $(v - r)^4/\eta_0^4$. Mello, Maciel and Zanchin, Phys. Rev. D 95, 084031, describe the conformal factor of this spacetime as a function of an advanced time of that kind.

Faraoni's metric in the cosmic time and the areal radius, in his section III, is not a chart of this spacetime: it follows from his (2.3), whose first line takes the scale factor for a function of Schwarzschild's time alone.
A footnote of Carrera and Giulini's section IV finds the step, and Mello, Maciel and Zanchin identify the metric it leads to as the nonrotating Thakurta metric, whose Ricci scalar does diverge on $r = r_s$.
For Sultana and Dyer's metric $R = (12\eta_0^4/\eta^6)(1 + r_s/r - r_s\eta/r^2)$, finite on the horizon.

## Step 3: what is checked

`sultana_dyer_check` holds the Kerr-Schild chart to:

- $a^2$ times the published metric of `schwarzschild`, pulled back along $ct = \eta - r_s\ln(r/r_s - 1)$;
- the two fluids as Saida, Harada and Maeda write them, at $M = r_s/2$, $G_{\mu\nu} = \mu\,u_\mu u_\nu + \tau\,k_\mu k_\nu$ with $u$ a unit timelike vector, $k$ null and $k\cdot u = -1$, $\mu = 12\eta_0^4D/(r^2\eta^6)$ and $\tau = \eta_0^4r_s(8r^2 + 3r_s(2r - \eta))/(r^2\eta^5D)$, $D = r^2 + r_s(r - \eta)$;
- the published conformal chart of `frw` at $k = 0$ and $a = \eta^2/\eta_0^2$, when $r_s = 0$;
- the Ricci scalar above;
- $|\nabla R|^2 = 0$ for the areal radius $R = \eta^2r/\eta_0^2$ on $r = \eta/2$ and on $r = (\sqrt{\eta^2 + 12r_s\eta + 4r_s^2} - \eta - 2r_s)/4$, Saida, Harada and Maeda's two trapping horizons, their $r_2$ and $r_1$.

Each other chart is held to being the Kerr-Schild chart pulled back.
`_tools/test_sultana_dyer.py` holds the published components to the same facts without sympy's tensors.

The sign of $D$ is the sign of the dust's density, so the energy conditions hold where $\eta < r(r + r_s)/r_s$, as Saida, Harada and Maeda and Majhi state; on the horizon that is $\eta < 2r_s$.

## Step 4: what the diagrams draw

Every diagram takes $r_s = 1$ and $\eta_0 = 3\,r_s$.
The null rays and the conformal diagram do not depend on $\eta_0$; the embedding does, through $a$.

- Spacetime diagrams: one plane for each chart. The Kerr-Schild plane marks $g^{rr} = 0$, the event horizon, and both trapping horizons. The plane of Schwarzschild's time declares the big bang as $(r/r_s - 1)e^{ct/r_s} = 1$, which is $\eta = 0$ in a form finite on $r_s$, so that the curve is found up to the top of the drawing. The plane of $v$ and $r$ sets `null_radius`, since $r$ at constant $v$ is a null direction, and `other_trapped`, which marks the trapped spheres of the outgoing family as well. `--verify` checks the rays against $\eta + r$ and $\eta - r - 2r_s\ln|r/r_s - 1|$, $ct \pm r_*$, and $v$ and $v - 2r_*$.
- Conformal diagram: Kruskal and Szekeres's map of Schwarzschild's ingoing chart, $V = e^{(\eta + r)/2r_s}$ and $U = (1 - r/r_s)e^{(r - \eta)/2r_s}$, restricted to $\eta > 0$. The big bang is the curve from $(U, V) = (1, 1)$, on the singularity, to $i^0$, so the region is a triangle with corners $(X, T) = (0, \pi/2)$, $i^+ = (\pi/2, \pi/2)$ and $i^0 = (\pi, 0)$, the shape Saida, Harada and Maeda's figure 3 gives its regions I and II. One view for each chart.
- Embedding diagram: the equator of a moment of $\eta$ in the Kerr-Schild chart, $a^2((1 + r_s/r)\,dr^2 + r^2d\phi^2)$, which is the paraboloid $z^2 = 4ar_s\rho$ over the areal radius $\rho = ar$, the Kerr-Schild slice of a Schwarzschild hole of radius $ar_s$. The geometry evolves, so it is a movie, $\eta$ from $3$ to $6\,r_s$, each frame out to $r = 4\,r_s$ and each the first one magnified by $a$.

## Step 5: the names in the History

Each name is as its paper's author line prints it, in Crossref's record of the paper.
Megan McClure, whose paper prints M L McClure, is as INSPIRE's author record 1047240 gives her, at the University of Toronto's astronomy department, Dyer's.
Thakurta's paper, Indian J. Phys. 55B, 304 (1981), has no record in Crossref or INSPIRE, so the History cites Mello, Maciel and Zanchin for his metric and gives no first name.
