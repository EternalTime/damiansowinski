# Hayward

Sean Hayward's regular black hole is his Eqs. (1) and (5), in his own two lengths,

$$ds^2 = -F\,c^2dt^2 + \frac{dr^2}{F} + r^2\,d\Omega^2,\qquad F = 1 - \frac{2mr^2}{r^3 + 2m\ell^2},$$

with $m = GM/c^2$ the mass and $\ell$, his $l$, the length of the core.
Its four charts are written by `_tools/derivations/print_charts.py --metric hayward`, and `verify_metrics.py --system hayward/<chart>` checks each in seconds.
Equation numbers are those of arXiv:gr-qc/0506126v2, the text of Physical Review Letters 96, 031103.

## Step 1. The parameters and the charts

The parameters are Hayward's, and no factor of $G$ or $c$ stands in a component: $m$ and $\ell$ are lengths, as Majumdar and Papapetrou's $m$ and Taub-NUT's $m$ and $l$ are.
Schwarzschild's $r_s$ is $2m$, and the metric is Schwarzschild's at $\ell = 0$ and flat at $m = 0$, his remark under (5).
Each chart has a source in the paper:

- `static`, his (1) with (5): the static time $t$ and the area radius $r$.
- `eddington_finkelstein_ingoing`, his (10) and (11): the advanced time $v = t + \int dr/F$, in which $ds^2 = r^2dS^2 + 2\,dv\,dr - F\,dv^2$.
- `eddington_finkelstein_outgoing`, his (20) at constant mass: the retarded time $u$, in which $ds^2 = r^2dS^2 - 2\,du\,dr - F\,du^2$.
- `evaporating`, his (11) with "the mass to depend on advanced time, $m(v)$, defining $F(r, v)$ by the same expression (5)": the chart his black hole forms and evaporates in.

He calls none of them by Eddington's or Finkelstein's name; the ids follow the charts of the same construction the collection already has for Schwarzschild, Kottler, and the rest.
A Kruskal chart gives $r$ only implicitly, so none is printed; the conformal diagram draws its own map through both horizons.

## Step 2. The centre, infinity, and the Einstein tensor

`hayward_check` in `print_charts.py` holds every chart to three things before it is written.

At the centre $F = 1 - r^2/\ell^2 + r^5/(2m\ell^4) + \dots$, his (3): de Sitter's static metric with $\Lambda = 3/\ell^2$, his (4), for every mass.
The term in $r^5$ is odd, so the metric is not analytic in $r^2$ at the centre.
Far away $F = 1 - 2m/r + O(r^{-4})$, his (2).

The mixed Einstein tensor is his (7), (8), and (12),

$$G^t{}_t = G^r{}_r = -\frac{12\ell^2m^2}{(r^3 + 2\ell^2m)^2},\qquad G^\theta{}_\theta = G^\phi{}_\phi = \frac{24(r^3 - \ell^2m)\ell^2m^2}{(r^3 + 2\ell^2m)^3},\qquad G^r{}_v = \frac{2r^4m'}{(r^3 + 2\ell^2m)^2},$$

the last only where the mass is a function of $v$, and with the opposite sign as $G^r{}_u$ in retarded time.
With $8\pi T = G$ the density is $\rho = -G^t{}_t/8\pi \ge 0$, the radial pressure is $-\rho$, and the transverse pressure is $p_\perp = \rho\,(2r^3 - 2\ell^2m)/(r^3 + 2\ell^2m)$.
So $\rho + p_\perp = 3\rho r^3/(r^3 + 2\ell^2m) \ge 0$ and the weak energy condition holds, $\rho + p_r + 2p_\perp = 2p_\perp$ is negative inside $r^3 = \ell^2m$ and the strong one fails there, and $p_\perp > \rho$ beyond $r^3 = 4\ell^2m$, where the dominant one fails.
The Weyl tensor carries the factor $r^3(r^3 - 4m\ell^2)$, so it vanishes at the centre and on that same sphere.
At the centre the Ricci scalar is $12/\ell^2$ and the Kretschmann scalar $24/\ell^4$, de Sitter's.

## Step 3. The horizons

$(r^3 + 2m\ell^2)F = r^3 - 2mr^2 + 2m\ell^2$, a cubic with one negative root and, for $\ell < \ell_* = 4m/3\sqrt{3}$, two positive ones, $r_- < r_+$.
At $\ell_*$ they meet at $r = 4m/3 = \sqrt{3}\,\ell_*$, which is his critical mass $m_* = 3\sqrt{3}\,\ell/4$ and radius $r_* = \sqrt{3}\,\ell$.
Writing the roots as $2m$ times $a$, $b$, and $1 - a - b$, the vanishing of the term in $r$ gives $ab = s(s - 1)$ with $s = a + b$, so the roots are rational when $s(4 - 3s)$ is the square of a rational.
Short of the extremal $s = 4/3$, the one with the smallest denominator is $s = 9/7$, with $a, b = 3/7, 6/7$:

$$\ell = \frac{12m}{7\sqrt{7}} = 0.648\,m,\qquad r^3 - 2mr^2 + 2m\ell^2 = \left(r - \tfrac{12m}{7}\right)\left(r - \tfrac{6m}{7}\right)\left(r + \tfrac{4m}{7}\right).$$

Every diagram is drawn there: $r_+ = 12m/7$, $r_- = 6m/7$, surface gravities $\kappa_+ = F'(r_+)/2 = 1/6m$ and $\kappa_- = 5/12m$, and $\ell/\ell_* = 0.84$.

## Step 4. The tortoise coordinate

$1/F = 1 + 2mr^2/\left((r - r_+)(r - r_-)(r - r_0)\right)$ is $1$ plus three simple poles with residues $A_i = 1/F'(r_i)$, so

$$r_* = r + \sum_i A_i\ln\left|1 - \frac{r}{r_i}\right| = r + 3m\ln\left|1 - \frac{7r}{12m}\right| - \frac{6m}{5}\ln\left|1 - \frac{7r}{6m}\right| + \frac{m}{5}\ln\left(1 + \frac{7r}{4m}\right)$$

vanishes at $r = 0$, which is how the charts' conventions fix its constant.
The negative root is a root of the cubic and no horizon, and the conformal diagram's `Tower` takes it with the other two, since the three residues together make $1/F - 1$.
`null_rays.py --verify` holds the traced rays of the static and both Eddington-Finkelstein planes to $ct \pm r_*$.

## Step 5. Printing

`chart_printer.Printer` writes a product as the line element does, the mass and the length before the radius, $2m\ell^2$ and $2mr^2$, and a sum by falling powers of $r$, so every value is printed around $r^3 + 2m\ell^2$ and $r^3 - 2mr^2 + 2m\ell^2$, the denominator and the numerator of $F$.
In the forming and evaporating chart each value is collected by $\partial_vm$, which sets what the radiation adds apart from what the static chart has already, as in $\Gamma^r{}_{vv}$ and $G_{vv}$.

## Step 6. The mass that forms and evaporates

Hayward asks for a profile with $m'$ at least continuous that is zero, rises, holds at $m_0 > m_*$, falls, and is zero again, his (13) to (17).
Every diagram of the `evaporating` chart takes, in units of $m_0$,

$$m(v) = \sin^2\frac{\pi v}{4}\ \ (0 \le v \le 2),\qquad 1\ \ (2 \le v \le 4),\qquad \cos^2\frac{\pi(v - 4)}{8}\ \ (4 \le v \le 8),$$

and zero outside, `HAYWARD_MASS` in `null_rays.py`, with the $\ell$ of Step 3.
Its derivative is continuous and vanishes at each join.
Trapped spheres exist while $m > m_* = 0.842\,m_0$, from $v = (4/\pi)\arcsin\sqrt{m_*} = 1.479$ to $v = 4 + (8/\pi)\arccos\sqrt{m_*} = 5.042$, his $v_b$ and $v_e$ of (18) and (19), and the curve $g^{rr} = 0$ opens and closes at $r = \sqrt{3}\,\ell = 1.122\,m_0$.
The drawings model the ingoing radiation alone: his outgoing radiation beyond a pair creation surface $r_0$, (20) to (23), changes nothing inside $r_0$, where the trapped region lies.

## Step 7. The conformal diagrams

The static hole is a `Tower` of the three roots, Reissner-Nordström's with $r = 0$ a regular centre: $r_*(0) = 0$ puts it on $u = -v$, the vertical lines $X = \pm\pi/2$, his Fig. 1.
`hayward` in `conformal.py` draws it symmetric about one pair of exteriors, a cell $k$ periods up the tower being the cell's own map moved by $k\pi$ in both $p$ and $q$.
The static chart covers the exterior I, the black hole II, and the inner region III.
The ingoing chart covers I, II, and the inner region III′ on the other side, with the cell's own time $t = r_* - v$ there, since a ray of constant $v$ keeps $q = \arctan e^{\kappa_+v}$ through all three.
The outgoing chart covers I, the white hole IV, and the inner region III′ one period down, with $t = -u - r_*$ there.
Across $r_-$ the map is continuous and no more, each ray's own coordinate closing on its limit as $(r - r_-)^{\kappa_+/\kappa_-}$, the power $2/5$.

The hole that forms and evaporates has the causal structure of flat space, his Fig. 5 inside $r_0$.
Every outgoing ray leaves the regular centre, at an advanced time $v_0$, so $p = \Phi(v_0)$ and $q = \Phi(v)$ with $\Phi(w) = \arctan((w - 4)/4)$ put the centre on $X = 0$ and the spacetime on Minkowski's triangle.
`hayward_left_centre` finds $v_0$ by integrating $dr/dv = F/2$ backward with the classical Runge-Kutta rule in steps of $0.004$ that end on the joins of the mass, and finishing each ray on the zero of its last step's Hermite cubic.
The script checks the map against the published metric at eight thousand events, sends six rays out again from the centre to within $10^{-7}$ of the events they were traced from, and follows forty-one rays from the centre to $r > 90\,m_0$, which is the absence of an event horizon.

## Step 8. The embedding

On the equator of $t = 0$ the metric is $dr^2/F + r^2d\phi^2$, and $F \le 1$ wherever it is positive, so flat space carries the whole slice at $dz/dr = \sqrt{1/F - 1}$.
Outside $r_+$ it is two sheets through the outer bifurcation sphere, as Schwarzschild's.
Inside $r_-$ it leaves the centre as the sphere of radius $\ell$ does, $z = r^2/2\ell + O(r^4)$, stands vertical at $r_-$, and runs on into the second inner region: a closed surface.

A slice of constant $v$ is null, so the forming and evaporating hole is drawn on slices of constant $v - r = T$, as Vaidya's shell is, on which the metric pulls back to $(2 - F)\,dr^2 + r^2d\phi^2$ with $F$ at $v = T + r$, and $dz/dr = \sqrt{2mr^2/(r^3 + 2m\ell^2)}$.
Each frame is checked against an independent quadrature of that slope, and the circles where $g^{rr}$ vanishes on it are found by bisection.
