# The universe on a wall moving past a black hole in five dimensions

The charts of `brane_wall_universe.json` are written by `print_charts.py --metric brane_wall_universe`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, the closed form of the Gaussian normal chart, the junction condition that moves the wall, and the wall the diagrams draw.

## Step 1. The spacetime

A wall with three dimensions of space in anti-de Sitter space of five, with the mirror symmetry of Randall and Sundrum's wall (`randall1999b`): two copies of one region glued along the wall.
With the symmetries of a Friedmann universe on the wall, the region is a piece of the black hole of anti-de Sitter space in five dimensions,

$$ds^2 = -f\,c^2dt^2 + \frac{dr^2}{f} + r^2d\Sigma_k^2, \qquad f = k + \frac{r^2}{\ell^2} - \frac{\mu}{r^2},$$

Kraus's (8) (`kraus1999`), with $d\Sigma_k^2 = d\chi^2/(1 - k\chi^2) + \chi^2d\Omega^2$ as Bowcock, Charmousis and Gregory write the unit 3-space (`bowcock2000`, their (2)).
Bowcock, Charmousis and Gregory prove that no other region borders such a wall when the bulk has a cosmological constant.
The wall is $r = a(\tau)$, and each copy is kept where $r$ falls away from the wall, so each holds the horizon $r_h$, the largest root of $f$, $r_h^2 = (\ell^2/2)(\sqrt{k^2 + 4\mu/\ell^2} - k)$ (Kraus's (9)).

- `static`: Kraus's chart, $r_h < r \le a(t)$, with the wall's radius declared as a function of the static time.
- `gaussian_normal`: the chart of Binétruy, Deffayet, Ellwanger and Langlois (`binetruy2000brane`, their (2) with $\dot b = 0$, (24) and (26)), in which the wall stands at $w = 0$ and $w$ is the proper distance from it, in the closed form of Mukohyama, Shiromizu and Maeda (`mukohyama2000global`, their (23) and (26), lower signs).

`brane_wall_universe_check` holds both charts to $R_{\mu\nu} = -4g_{\mu\nu}/\ell^2$ and the static chart to $K = 40/\ell^4 + 72\mu^2/r^8$.

## Step 2. Why each chart covers one side

The two copies are glued with $\partial_w g_{ij}$ changing sign, so a chart that ran through the wall would carry the wall's matter as a delta in its curvature.
In the Gaussian normal chart that is $|w|$ in place of $w$, and every name of the chart would hold $|w|$ and its derivatives, $\mathrm{sgn}(w)$ and $\delta(w)$.
Each chart is therefore one copy, ending on the wall, as the one sided Poincaré chart of `randall_sundrum` is, and the wall's matter enters through the junction condition of Step 4.
`OVERRULED` in `metric_tags.py` says so for the tags: each chart is an Einstein space, and the spacetime, with the wall's matter on $w = 0$, is not.

## Step 3. The Gaussian normal chart

Mukohyama, Shiromizu and Maeda follow the spacelike geodesics normal to the wall out of the static chart.
With $E = \dot a$ conserved along each, $u^\mu\partial_\mu = (E/f)\partial_T - \sqrt{f + E^2}\,\partial_r$ on the side where $r$ falls, and

$$2r^2 + \ell^2(E^2 + k) = \sqrt{4\ell^2\mu + \ell^4(E^2 + k)^2}\cosh\left(\frac{2(w_0 - w)}{\ell}\right),$$

their (7), with $w_0(t)$ set by $r = a$ at $w = 0$.
With $t$ the proper time on the wall and a dot for $d/d(ct)$, $r^2 = a^2\varphi$ with

$$\varphi = \cosh\frac{2w}{\ell} + \frac{\ell^2(\dot a^2 + k)}{2a^2}\left(\cosh\frac{2w}{\ell} - 1\right) - Q\sinh\frac{2w}{\ell}, \qquad Q = \sqrt{1 + \ell^2\left(\frac{\dot a^2 + k}{a^2} - \frac{\mu}{a^4}\right)},$$

and the metric is $-(\psi^2/\varphi)\,c^2dt^2 + dw^2 + a^2\varphi\,d\Sigma_k^2$ with $\psi = \varphi + \partial_t\varphi/(2H)$, $H = \dot a/a$.

Writing $X = \ell^2(\dot a^2 + k)/2a^2$, $\partial_tX/(2H) = \ell^2(a\ddot a - \dot a^2 - k)/2a^2$, so $X + \partial_tX/(2H) = \ell^2\ddot a/2a$.
From $Q^2 = 1 + 2X - \ell^2\mu/a^4$, $Q + \partial_tQ/(2H) = \left(2a^2 + \ell^2(\dot a^2 + a\ddot a + k)\right)/2a^2Q$.
So

$$\psi = \cosh\frac{2w}{\ell} + \frac{\ell^2\ddot a}{2a}\left(\cosh\frac{2w}{\ell} - 1\right) - \frac{2a^2 + \ell^2(\dot a^2 + a\ddot a + k)}{2a^2Q}\sinh\frac{2w}{\ell},$$

which is their (26) with $k_{\text{MSM}} = 2/\ell$, $1 + 2k_{\text{MSM}}^{-2}(2H^2 + \dot H + k/a^2)$.
The arXiv copy of their (24) prints the same factor as $1 + \ell^2(2H^2 + \dot H + k/a^2)$, without the half that (26) and the definition of $\psi$ carry.
The lapse $N = \psi/\sqrt\varphi$ is $\partial_tA/\dot a$ for the areal radius $A = a\sqrt\varphi$, which is Binétruy, Deffayet, Ellwanger and Langlois's (38).

`brane_wall_universe_check` holds every component of the chart's Riemann tensor, with $A$, $N$ and $Q$ written out, to the bulk's own form, the Weyl tensor of the black hole beside the curvature of anti-de Sitter space,

$$R_{abcd} = \frac{\mu}{A^4}\left(\tfrac{3}{2}P \owedge P - P \owedge S + \tfrac{1}{2}S \owedge S\right)_{abcd} - \frac{1}{\ell^2}\left(g_{ac}g_{bd} - g_{ad}g_{bc}\right),$$

with $P$ the metric of the plane of $t$ and $w$, $S$ that of the 3-space and $\owedge$ the Kulkarni-Nomizu product, and the Riemann, Ricci and Einstein tensors are printed in that form.
The Weyl tensor and the Kretschmann scalar are printed as the names have them: written out, one Weyl component took over two minutes to compare, and the checker compares them in the names at once.

## Step 4. The junction condition

With the mirror symmetry, $[K_{\mu\nu}] = 2K_{\mu\nu}$ on the side $w > 0$, and Israel's conditions (`israel1966`) for a wall of energy density $\rho_b = \sigma + \rho$ and pressure $p_b = p - \sigma$ are those of Binétruy, Deffayet, Ellwanger and Langlois (their (19) and (20)),

$$2K^i{}_i = \frac{2\,\partial_wA}{A} = -\frac{\kappa^2\rho_b}{3}, \qquad 2K^t{}_t = \frac{2\,\partial_wN}{N} = \frac{\kappa^2(2\rho_b + 3p_b)}{3},$$

at $w = 0$, with $\kappa^2 = 8\pi G_5/c^4$.
At the wall $A = a$, $N = 1$ and $\partial_wA = -aQ/\ell$, so the first is $Q = \kappa^2\ell\rho_b/6$, which is $Q = 1 + \rho/\sigma$ for Randall and Sundrum's tension $\kappa^2\sigma = 6/\ell$, $\sigma = 3c^4/4\pi G_5\ell$.
Squared, it is the Friedmann equation the History states,

$$H^2 = \frac{8\pi G_4}{3c^2}\rho\left(1 + \frac{\rho}{2\sigma}\right) - \frac{kc^2}{a^2} + \frac{\mu c^2}{a^4}, \qquad G_4 = \frac{G_5}{\ell},$$

with $H = (da/d\tau)/a$, Kraus's (26) with $\mu_+ = \mu_-$ (`kraus1999`) and Ida's (28) with $\alpha_R = \alpha_L$ (`ida2000`).
The second, with the first, is $\dot\rho_b + 3H(\rho_b + p_b) = 0$, which is the conservation of the matter on the wall.
In the static chart the same condition is $\sqrt{f(a) + \dot a^2} = (a/\ell)(1 + \rho/\sigma)$, the energy of the wall along the static Killing vector, so $f(a)\,dt/d\tau = (a/\ell)(1 + \rho/\sigma)$.
`_tools/test_brane_wall_universe.py` reads both conditions off the published Gaussian normal chart and holds them to the Friedmann equation and to conservation.

## Step 5. The wall drawn

Every diagram draws the wall with no matter beside its tension, $\rho = 0$, in a flat space, $k = 0$, past a black hole of $\mu = \ell^2$, in units of $\ell$, Kraus's unbounded motion at $\sigma = \sigma_c$, $k = 0$.
Then $f = r^2 - 1/r^2$, the horizon is $r_h = 1$ with surface gravity $f'(1)/2 = 2$, and $\dot a^2 = 1/a^2$, so $a^2 = 2c\tau$, the expansion of a universe of radiation, all of it the black hole's.
With the tortoise radius $r_* = \tfrac14\ln|(r - 1)/(r + 1)| + \tfrac12\arctan r$, $dt/da = a^4/(a^4 - 1)$ integrates to

$$ct = a + r_*(a) - \arctan a, \qquad u = a - \arctan a, \qquad v = a + 2r_*(a) - \arctan a$$

on the wall, so in Kruskal's $U = -e^{-2u}$, $V = e^{2v}$, with $UV = (1 - r)e^{2\arctan r}/(1 + r)$, the wall is $U = -e^{-2(a - \arctan a)}$ and $V = -(1 - a)e^{2a}/(1 + a)$ for every $a > 0$.
It leaves $r = 0$ at $U = V = -1$, the middle of the past singularity, crosses the past horizon $V = 0$ at $U = -e^{\pi/2 - 2}$, and reaches the corner of the boundary and the future singularity as $a \to \infty$, staying below the future horizon $U = 0$.
In the Gaussian normal chart $Q = 1$, and

$$A^2 = a^2\varphi, \qquad \varphi = e^{-2w} + \frac{\sinh^2 w}{4c^2t^2}, \qquad N = \frac{\psi}{\sqrt\varphi}, \qquad \psi = e^{-2w} - \frac{\sinh^2 w}{4c^2t^2},$$

so the lapse vanishes on $w_* = \tfrac12\ln(1 + 4ct)$.
`brane_wall_universe.py` holds these forms.

## Step 6. The diagrams

- The spacetime diagram of the static chart is its plane of $t$ and $r$ at $\chi = 1/2$, $\theta = \pi/2$, $\phi = 0$, from $r = 0$ to $3\ell$, the region $r > a(t)$ beyond the wall hatched by `where`, with the wall drawn by `curves` from the declared $a(t)$, the root of the closed form of Step 5. `--verify` holds its rays to $ct \pm r_*$.
- The spacetime diagram of the Gaussian normal chart is its plane of $t$ and $w$ for $0 < ct < 2\ell$ and $0 \le w \le 2\ell$, with $a = \sqrt{2ct\ell}$ declared, hatched beyond $w_*$, which is drawn.
- The conformal diagram is one side of the wall in the planar black hole's Kruskal square, by $p = \arctan U$ and $q = \arctan V$ of Step 5: the boundaries lie on $\tan p\tan q = -e^\pi$ and the singularities on $T = \pm\pi/2$. The static chart's four cells are checked against the published metric, the wall against its static time, and the region beyond the wall is drawn as removed.
- The embedding diagram is the wall's own moment $t = \ell/c$, the surface of $w$ and $\chi$, rolled so that $\chi$ runs once round, $\rho = A$. Along it $2r^2 + E^2 = \sqrt{4 + E^4}\cosh(2(w_0 - w))$ with $E = 1/a$ (Mukohyama, Shiromizu and Maeda's (7)), so the moment falls from the wall at $r = \sqrt2$ through the horizon to $r_{\min} = 0.884$ at $w_0 = 0.708$, inside the black hole, and grows through the other horizon into the other exterior. $(\partial_wA)^2 = A^2 - 1/A^2 + 1/2$ exceeds 1 outside the two circles $r = 1.1317$, where the surface is drawn in Minkowski space.
- The moment is drawn on the other diagrams through `brane_wall_universe.Moment`: on the Gaussian normal plane as the line $t = \ell/c$ out to $w_*$; on the static plane from the wall to the horizon and through the black hole from horizon to horizon, with $t = (\ln|V| - \ln|U|)/4$; and on the conformal diagram by $U$ and $V$, followed from the wall by $dv/dw = -1/(E + \sqrt{f + E^2})$ to $w_0$ and by $du/dw$ of the same form beyond it, both smooth through the horizons.
