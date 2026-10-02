# The Boulware-Deser black hole

Boulware and Deser's solution of Einstein-Gauss-Bonnet gravity in five dimensions is $ds^2 = -f\,c^2dt^2 + dr^2/f + r^2d\Omega_3^2$ with

$$f_\mp = 1 + \frac{r^2}{2\ell^2}\left(1 \mp \sqrt{1 + \frac{4\ell^2r_0^2}{r^4}}\right).$$

Its four charts are written by `_tools/derivations/print_charts.py --metric boulware_deser`, and `verify_metrics.py --system boulware_deser/<chart>` checks each in seconds.

## Step 1: the parameters

Torii and Maeda's action is $(R + \alpha L_{GB})/2\kappa_n^2$ with $L_{GB} = R^2 - 4R_{\mu\nu}R^{\mu\nu} + R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$, and their solution, (9) of Phys. Rev. D 71, 124002, is $f = 1 + (r^2/2\tilde\alpha)(1 \mp \sqrt{1 + 4\tilde\alpha\tilde M/r^{n-1}})$ with $\tilde\alpha = (n-3)(n-4)\alpha$ and $\tilde M = 16\pi GM/((n-2)\Sigma_{n-2})$.
At $n = 5$ both carry the dimension of an area, so the charts quote two lengths: $\ell^2 = \tilde\alpha = 2\alpha$ and $r_0^2 = \tilde M = 8GM/(3\pi c^2)$, which is the square of the horizon radius of Tangherlini's black hole of the same mass.
Garraffo and Giribet's (21), Mod. Phys. Lett. A 23, 1801, is the same function with $4\alpha$ for $2\ell^2$ and $2M$ for $r_0^2$.

## Step 2: the radical

The charts name $W = \sqrt{r^4 + 4\ell^2r_0^2}$, an area.
Since $W^2 - r^4 = 4\ell^2r_0^2$, $f_- = 1 - 2r_0^2/(r^2 + W)$ and $f_+ = 1 + 2r_0^2/(W - r^2)$, so the second branch is the first with $W \to -W$.
`boulware_deser_pretty` writes every value in $r$, $r_0$ and $W$, with $\ell^2 = (W^2 - r^4)/4r_0^2$, where no relation is left among the generators.
In that form each value of the first branch is Tangherlini's at $W = r^2$, which is $\ell = 0$: the Ricci tensor carries the factor $W - r^2$ and the Kretschmann scalar becomes $72r_0^4/r^8$.

## Step 3: the sources of the charts

- `spherical`: Boulware and Deser's own, Phys. Rev. Lett. 55, 2656, the branch that is asymptotically Schwarzschild.
- `eddington_finkelstein_ingoing` and `eddington_finkelstein_outgoing`: $v, u = ct \pm r_*$ with $dr_*/dr = 1/f$, Kobayashi's (12), Gen. Rel. Grav. 37, 1869, with $\epsilon = \pm 1$, and Charmousis's (54), Lect. Notes Phys. 769, 299.
- `spherical_plus`: the other root, which Boulware and Deser exhibit beside the first and which approaches anti-de Sitter space of radius $\ell$.

`boulware_deser_check` holds every chart to $G_{\mu\nu} + (\ell^2/2)H_{\mu\nu} = 0$ with Torii and Maeda's $H_{\mu\nu}$, their (5), at three random points to forty digits; the first branch to its horizon $r_h^2 = r_0^2 - \ell^2$, to the surface gravity $r_h/(r_h^2 + 2\ell^2)$, to $f(0) = 1 - r_0/\ell$ and to Tangherlini's metric as $\ell \to 0$; each Eddington-Finkelstein chart to being the static one pulled back; and both branches to $r^4K \to 12r_0^2/\ell^2$ at the centre, Torii and Maeda's (28).
`_tools/test_boulware_deser.py` runs the same check on the metric components in the file.

## Step 4: the tortoise coordinate

With $W_h = r_h^2 + 2\ell^2$ the value of $W$ on the horizon,

$$\frac{1}{f} = \frac{r^2 + 2\ell^2 + W}{2(r^2 - r_h^2)} = 1 + \frac{W_h}{r^2 - r_h^2} - \frac{W - r^2 + 2\ell^2}{2(W + W_h)},$$

so $r_* = r + (W_h/2r_h)\ln|(r - r_h)/(r + r_h)| - \tfrac12 J(r)$, Tangherlini's with the surface gravity $r_h/W_h$ in place of $1/r_h$, less half the integral $J$ from $0$ of $(W - s^2 + 2\ell^2)/(W + W_h)$, which is smooth, falls as $2\ell^2/s^2$ and vanishes at $\ell = 0$.
`slices.boulware_deser_rstar` sums $J$ by Gauss and Legendre's rule, in $s$ to $s = 1$ and in $1/s$ beyond.
The other branch has $1/f_+ = 2\ell^2/(r^2 + 2\ell^2 + W)$, so its $r_*$ runs from $0$ at the centre to a finite $R$ at infinity, `slices.boulware_deser_plus_rstar`.

## Step 5: what the diagrams draw

The black hole is drawn at $r_0 = 13r_h/12$ and $\ell = 5r_h/12$, in units of its horizon radius, where $\kappa = 72/(97\,r_h)$ and $f(0) = -8/5$.
Beroiz, Dotti and Gleiser, Phys. Rev. D 76, 024012, find the hole unstable for $3/2 < \mu/\alpha < 9/2 + 3\sqrt2$, their (31), which with $r_0^2 = 2\mu/3$ and $\ell^2 = \alpha$ is $1 < r_0/\ell < 1 + \sqrt2$; the ratio drawn, $13/5$, lies above that range.
The other branch is drawn at $r_0 = \ell$, where $f_+(0) = 2$ and $R = 1.1981\,\ell$.

- Spacetime diagrams: every chart, the Eddington-Finkelstein charts also against their own null coordinate; `--verify` checks the rays against $ct \pm r_*$.
- Conformal diagrams: the black hole is Kruskal and Szekeres's square by `BoulwareDeserTower`, with $r_*(0) = 0$ so that the singularity lies on $T = \pm\pi/2$, spacelike, as Torii and Maeda find for $\tilde M > \tilde\alpha$; the other branch is the strip $X = \pi r_*/2R$, $T = \pi ct/2R$ between a timelike singularity and a timelike boundary.
- Embedding diagrams: the static slice of the black hole in flat space, $dz/dr = \sqrt{2r_0^2/(r^2 + W - 2r_0^2)}$, through the throat into the other exterior, $2.77\,r_h$ high at $6\,r_h$; and the other branch in Minkowski space, $dZ/dr = \sqrt{1 - 1/f_+}$, a sheet that leaves the centre along a cone of slope $\sqrt{r_0/(r_0 + \ell)}$. The slices are static, so neither is a movie.

The branch's cone sets its parameter: the check of chords across a surface compares a straight line in space with a line that turns as it runs out, which near an apex differ by $f_+(0)\,\Delta\phi^2/6$ of their length at any step, so the check passes only for $f_+(0) < 3$.

## Step 6: the names in the History

Each name is as its paper's author line prints it; Stanley Deser, whose papers print S. Deser, and David Wiltshire, whose paper prints D. L. Wiltshire, are as INSPIRE's author records 1011955 and 983425 give them, and Lanczos is Kornel, as everywhere on the site.
