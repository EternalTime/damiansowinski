# Brill's charged Taub-NUT

Taub-NUT with a charge, Dieter Brill's solution of the Einstein-Maxwell equations of 1964,

$$ds^2 = -\frac{\Delta}{\Sigma}\left(c\,dt + 2l\cos\theta\,d\phi\right)^2 + \frac{\Sigma}{\Delta}dr^2 + \Sigma\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),$$

$$\Sigma = r^2 + l^2,\qquad \Delta = r^2 - 2mr - l^2 + r_q^2.$$

Its five charts are written by `_tools/derivations/print_charts.py --metric brill_charged_taub_nut`, about ten seconds each on 2 October 2026, and `verify_metrics.py --system brill_charged_taub_nut/<chart>` checks each in a few seconds.

## Step 1. Where the metric comes from

Brill's paper (`brill1964`) is behind the publisher's login, so it was read through its abstract on the publisher's page, its record at Crossref and the papers that restate it.
The abstract gives the solution as a closed universe of topology $S^3 \times R$, homogeneous and not isotropic, a generalisation of one of Taub's solutions.
Clément and Guenouche (`clement2018`) say it "was first given (in the Taub form) for $C = 0$ by Brill", and the form above is equations (2.1) and (2.2) of Clément, Gal'tsov and Guenouche (`clement2016`) and equations (1) and (2) of Halla and Perlick (`halla2023`), each at $C = 0$.
Their twist is written $-2n\cos\theta$ or $-2l\cos\theta$; the Taub-NUT page writes $+2l\cos\theta$, which is the same metric with $\phi$ reversed, and this page follows it, so that $r_q = 0$ is `taub_nut` term by term.
The mass enters as the length $m = GM/c^2$ and the charge as the charge radius $r_q$ of the Reissner-Nordström page, so that $l = 0$ is `rn_metric` at $r_s = 2m$.
The metric depends on the electric and magnetic charges only through the sum of their squares (`clement2016`).

## Step 2. The checks

`brill_charged_taub_nut_check` holds the first chart to the Einstein-Maxwell equations with the potential of a purely electric charge,

$$A = \frac{r_q\,r}{\Sigma}\left(c\,dt + 2l\cos\theta\,d\phi\right),$$

`clement2016`'s (2.2) at $p = 0$: $G_{\mu\nu} = 2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac{1}{4}g_{\mu\nu}F^2\right)$ in units where the charge is the length $r_q$, and $\partial_\mu\left(\sqrt{-g}\,F^{\mu\nu}\right) = 0$ with $\sqrt{-g} = \Sigma\sin\theta$.
It holds the same chart to the published metric of `taub_nut` at $r_q = 0$, of `rn_metric` at $l = 0$, and of the spherical chart of `israel_wilson_perjes` at $r_q^2 = m^2 + l^2$, where $\Delta = (r - m)^2$.
Every other chart is held to $R_{\mu\nu}R^{\mu\nu} = 4r_q^4/\Sigma^4$ and to being the first pulled back through `brill_charged_taub_nut_jacobian`.

The Ricci scalar vanishes, and the mixed Ricci tensor has the eigenvalues $\mp r_q^2/\Sigma^2$, each twice.
The Kretschmann scalar is written as `clement2016`'s (2.6), the square of the Ricci tensor plus the square of the Weyl tensor,

$$K = \frac{8r_q^4}{\Sigma^4} + \frac{48}{\Sigma^6}\Big(\left(m^2 - l^2\right)\left(r^6 - 15l^2r^4 + 15l^4r^2 - l^6\right) - 2mr\left(\left(r_q^2 - 6l^2\right)r^4 - 10l^2\left(r_q^2 - 2l^2\right)r^2 + l^4\left(5r_q^2 - 6l^2\right)\right) + r_q^2\left(\left(r_q^2 - 10l^2\right)r^4 - 2l^2\left(3r_q^2 - 10l^2\right)r^2 + l^4\left(r_q^2 - 2l^2\right)\right)\Big),$$

which the script checks against sympy and which is finite at every $r$ while $l \neq 0$.

## Step 3. How a value is written

Every chart declares $\Sigma$ and $\Delta$ as names, or $\Sigma = \tau^2 + l^2$ and $U = l^2 - r_q^2 + 2m\tau - \tau^2$ in the chart of Brill's universe.
`BrillForms` factors a value with the even powers of the sine written in the cosine, writes each factor that is a multiple of a name's polynomial by the name, $(\cos\theta + 1)(\cos\theta - 1)$ as $-\sin^2\theta$, and every other sum in the form with the fewest terms among: as it stands; with the mass written out by $2mr = r^2 - l^2 + r_q^2 - \Delta$; with the charge written out by $r_q^2 = \Delta - r^2 + 2mr + l^2$; each of those with $r^2 = \Sigma - l^2$; each in the cosine, in the sine, and with $A\cos^2\theta + B$ written $(A + B)\cos^2\theta + B\sin^2\theta$.
So $g_{\phi\phi} = \left(\Sigma^2\sin^2\theta - 4l^2\Delta\cos^2\theta\right)/\Sigma$ and $\Gamma^t{}_{tr} = \left(r\Sigma - m\Sigma - r\Delta\right)/(\Sigma\Delta)$, which is $\tfrac{1}{2}\partial_r\ln(\Delta/\Sigma)$.
The factoring names its generators in one order, the radius, $\cos\theta$, $\sin\theta$, $m$, $l$ and $r_q$.
Left to choose, sympy orders them by the run's hashes, and under the hash seed 26 one lowered Weyl component of Brill's universe did not factor in twenty minutes, where under the fixed order every chart prints in about ten seconds under that seed and others; a factor is recognised as a multiple of a name by one coefficient and one expansion, with no division.
The Weyl tensor's components are all multiples of $\Sigma^2 - \Sigma\Delta - r_q^2\Sigma + 4l^2\Delta$ and of $r\Sigma - m\Sigma - 2r\Delta$, the real and imaginary parts of the one Weyl scalar of a spacetime of Petrov type D.

## Step 4. The other four charts

The chart with one Misner string keeps $r$, $\theta$ and $\phi$ and takes $t_N = t + 2l\phi/c$, so that the twist becomes $-2l(1 - \cos\theta)$, which vanishes on the northern half of the axis.
It is the form of Newman, Tamburino and Unti's original metric (`nut1963`), Manko and Ruiz's $C = \mp 1$ as `halla2023` and `clement2016` write it, and the chart of the same name on the Kerr-Taub-NUT page at $a = 0$ with $r_q^2$ added to $\Delta$.

The Eddington-Finkelstein charts follow the two principal null congruences of that second form, $dv = c\,dt_N + \Sigma\,dr/\Delta$ and $du = c\,dt_N - \Sigma\,dr/\Delta$, with $r_* = 0$ at $r = 0$.
Then $g_{rr} = 0$ and each chart runs through both roots of $\Delta$ and through $r = 0$.
They are the Kerr charts of the Kerr-Taub-NUT page at $a = 0$, where no angle needs unwinding.

Brill's universe is the region between the horizons, $\tau_- < r < \tau_+$, in the time $\tau = r$ and the Euler angle $\psi = ct/2l$,

$$ds^2 = -\frac{\Sigma}{U}d\tau^2 + \frac{4l^2U}{\Sigma}\left(d\psi + \cos\theta\,d\phi\right)^2 + \Sigma\left(d\theta^2 + \sin^2\theta\,d\phi^2\right),\qquad U = -\Delta(\tau),$$

Hawking and Ellis's form of Taub's universe (`hawking1973`, section 5.8) with $r_q^2$ subtracted from their $U(t^2 + l^2)$, and the form Brill's abstract describes.
With $\psi$ of period $4\pi$ its moments are spheres $S^3$, which is Misner's periodic time, $ct \sim ct + 8\pi l$.

## Step 5. The values the diagrams are drawn at

The black hole is drawn at $m = 1$, $l = 3m/4$ and $r_q = m$: $\Delta = (r - 7/4)(r - 1/4)$, so the horizons are $r_+ = 7m/4$ and $r_- = m/4$, a charge that would make Reissner and Nordström's hole extremal with the twist keeping its horizons apart.
There $r_* = r + \tfrac{29}{12}\ln|1 - 4r/7| - \tfrac{5}{12}\ln|1 - 4r|$, which `slices.brill_rstar` sums.
The wormhole is drawn with no mass, in units of $l$, at $r_q = 3l/2$: $\Delta = r^2 + 5l^2/4$ has no root, the two sides of the throat are mirror images, and $r_* = r - \arctan(2r/\sqrt5)/(2\sqrt5)$.

The spacetime diagrams draw the equator's plane of $t$ and $r$ in the first chart, where $g_{t\phi}$ vanishes with $\cos\theta$ and the rays are the principal null rays; the regular half axis in the chart with one string and in both Eddington-Finkelstein charts; and the plane of $\tau$ and $2l\psi$ of Brill's universe.
The conformal diagram is the tower of the regular half axis, with the chart of one string, both Eddington-Finkelstein charts and Brill's universe, the cell between the horizons in which $r$ grows toward the future, on it, and the wormhole's diamond.
The equator is not totally geodesic when $l \neq 0$, so no conformal diagram of it is drawn; the embedding diagrams are its slices of constant $t$, the black hole's through the bifurcation sphere into a second exterior and the wormhole's through its throat, a surface of revolution of radius $\sqrt{r^2 + l^2}$ with $(dz/dr)^2 = \Sigma/\Delta - r^2/\Sigma$.
The wormhole's slice embeds at every $r$ because $r_q^2 \le 3l^2$ there: $\Sigma^2 - r^2\Delta = (3l^2 - r_q^2)r^2 + l^4$ at $m = 0$, which has to stay positive.
