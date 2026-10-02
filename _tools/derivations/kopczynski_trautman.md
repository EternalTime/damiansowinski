# Kopczyński and Trautman's universe with torsion

A flat universe of dust with every spin aligned, in the Einstein-Cartan theory, which turns round at a least size where Friedmann's dust has its big bang.
This note records the sources, the three charts, what each is checked against, and what the drawings draw.

## Step 1: the sources

W. Kopczyński, "A non-singular universe with torsion", Physics Letters A 39, 219 (1972), doi:10.1016/0375-9601(72)90714-1, and "An anisotropic universe with torsion", Physics Letters A 43, 63 (1973), doi:10.1016/0375-9601(73)90546-X; Andrzej Trautman, "Spin and Torsion May Avert Gravitational Singularities", Nature Physical Science 242, 7 (1973), doi:10.1038/physci242007a0.
Each record was checked on Crossref, and the two abstracts of Kopczyński were read on INSPIRE; his papers print his name as W. Kopczyński, and so does the page.
The full text of those three was not read: Nature's page gives the first paragraph and the reference list of Trautman's note, and the page quotes nothing else from it.
The line element and the field equation are taken from Trautman's own article "Einstein-Cartan Theory" in the Encyclopedia of Mathematical Physics (2006), arXiv:gr-qc/0606062, read in full, and the Friedmann equations with spin for any curvature from Nikodem Popławski, Physics Letters B 694, 181 (2010), arXiv:1007.0587, equations (10) to (15).

## Step 2: the scale factor

Trautman's section "Cosmology with spin and torsion" takes dust with $P^\mu = \rho u^\mu$, $u^\mu = \delta^\mu_0$ and the one spin component $S_{23} = \sigma$, the line element $dt^2 - R(t)^2(dx^2 + dy^2 + dz^2)$, and finds his equation (33) at $G = c = 1$,

$$\tfrac{1}{2}\dot R^2 - \frac{M}{R} + \frac{3}{2}\frac{S^2}{R^4} = 0, \qquad M = \tfrac{4}{3}\pi\rho R^3, \qquad S = \tfrac{4}{3}\pi\sigma R^3,$$

with $M$ and $S$ constant.
Multiplying by $18R^4$ gives $(d(R^3)/dt)^2 = 36MR^3 - 27S^2$, so $R^3 = 3S^2/2M + 9Mt^2/2$ with $t$ counted from the least radius.
With $a = R/R(0)$ and the units restored, $a^3 = 1 + c^2t^2/\ell^2$, where $\ell = c/\sqrt{6\pi G\rho_0}$ and $\rho_0$ is the density of the dust at $t = 0$.
The least radius is where $\rho c^4 = 2\pi G\sigma^2$, which is Popławski's effective energy density $\epsilon - \kappa s^2/4$ vanishing.
In these lengths $(da/d(ct))^2 = (4/9\ell^2)(1/a - 1/a^4)$.

Trautman's estimate, $10^{80}$ nucleons of mass $m$ with spins $\hbar/2$ aligned, gives $R(0)^3 = 3GN\hbar^2/8mc^4$, which is $R(0) = 1.27$ cm and a density of $1.9 \times 10^{55}$ g/cm$^3$, $10^{38.4}$ times below the Planck density; `_tools/test_kopczynski_trautman.py` holds the numbers the history states.

## Step 3: the charts

`comoving_cartesian` is Trautman's chart and `comoving_spherical` the same in spherical coordinates, which Kopczyński's first paper has by its abstract, "a spherically-symmetric gravitational field produced by spinning dust".
Both name $a = (1 + c^2t^2/\ell^2)^{1/3}$, and `print_charts.kopczynski_trautman` prints every value in $a$, $ct$ and $\ell$: the time is written $\ell s/c$, every fractional power of $1 + s^2$ becomes that power of $a^3$, every even power of $s$ a power of $a^3 - 1$, and each side of a value is factored in $a^3$.
`conformal` leaves $a(\eta)$ free, as the published conformal chart of `frw` does, since $\eta = \int c\,dt/a$ is a hypergeometric function of $t$ with no inverse in closed form; its parameter states the equation $a$ solves, $(\partial_\eta a)^2 = 4(a - a^{-2})/9\ell^2$.

`kopczynski_trautman_check` holds each chart to a vanishing Weyl tensor and to the Einstein tensor of the Christoffel connection that Trautman's effective stress tensor gives: $G^t{}_t = -(4/3\ell^2)(a^{-3} - a^{-6})$ and $G^x{}_x = -(4/3\ell^2)a^{-6}$, dust of density $\rho_0/a^3$ plus a stiff fluid of negative energy density.
The explicit scale factor is checked against the modified Friedmann equation, the spherical chart to being the Cartesian one pulled back, and the conformal chart on the equation its parameter states.
The published curvature is that of the metric alone; the torsion is not in it.

## Step 4: the drawings

The spacetime diagrams draw each chart through the turn in units of $\ell$.
The conformal chart's $a(\eta)$ is solved by `null_rays.DustSolver` from the chart's own published $G^r{}_r$ set to the pressure of the spins, $-4/3\ell^2a^6$, the solver's new `sources`, from $a = 1$ at rest; the rays of the comoving charts are checked against $\eta(t) \pm r$ with $\eta(t) = t\,{}_2F_1(1/3, 1/2; 3/2; -t^2)$.
The dotted curve is the Hubble sphere, $r = 3\ell^2a^2/2c|t|$, least at $ct = \pm\sqrt{3}\,\ell$, where it is $(\sqrt{3}/2)4^{2/3}\ell = 2.18\,\ell$.

The conformal diagram is the whole of Minkowski's half diamond, since $\eta$ grows as $3\ell(ct/\ell)^{1/3}$ without limit in both directions; Friedmann's flat dust fills only the upper half.
The embedding diagram is a movie of the equator of space from $ct = -3\ell$ to $3\ell$: every moment is a flat plane, and the disc of dust out to $r = \ell$ has the radius $a\ell$, least at the turn.
Its moments are marked on each spacetime diagram and on both conformal views, on the conformal chart at $\eta(t)$.
The acceleration $d^2a/d(ct)^2 = 2(4 - a^3)/9\ell^2a^5$ is positive until $a^3 = 4$, which is $ct = \sqrt{3}\,\ell$, the same moment the Hubble sphere comes nearest.
