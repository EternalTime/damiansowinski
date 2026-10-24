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
