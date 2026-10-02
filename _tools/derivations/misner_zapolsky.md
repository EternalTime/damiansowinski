# The star of infinite central density: the charts, the sources and what the diagrams draw

The spacetime is the static sphere of perfect fluid with $p = \rho c^2/3$ whose density is infinite at the centre, together with the two stars Tolman built on it and the family it belongs to.
This note records where each chart of `misner_zapolsky.json` comes from, what was read and what was not, and what the diagrams take from the charts.
Formulas here are in $G = c = 1$, as the sources write them.

## Step 1. What was read

Tolman (1939) and Oppenheimer and Volkoff (1939) were read in full, from the scans of the two papers.
Every equation number below is theirs.
Misner and Zapolsky's letter (1964) was verified on Crossref and INSPIRE (author line C. W. Misner and H. S. Zapolsky, *Physical Review Letters* **12**, 635-637) and its text was not available.
What the History says of it rests on Chavanis (2002), sections 3.6 and 4.1, who describes their model and their result, and on Dev and Gleiser (2002), where they take up the density that falls as $1/r^2$, who state the solution and its equation of state.
Chavanis (2002, 2008), Dev and Gleiser (2002), and the abstract of Sorkin, Wald, and Zhang (1981) were read from arXiv and INSPIRE.
Meltzer and Thorne (1966) is cited for the spiral of mass against radius on Chavanis's (2008) word, where he treats stars in $d$ dimensions and names their Fig. 2, and was not read.

Zapolsky's first name is from the Crossref record of his own paper of 1968 in *The Astrophysical Journal* (doi:10.1086/180244), which prints Harold S. Zapolsky, and Misner's from his INSPIRE author record.

Chavanis (2008) writes that the singular solution "was first found by Klein (1947)".
Tolman's section 8 and Oppenheimer and Volkoff's (22) hold it eight years earlier, so the History gives it to them.

## Step 2. The sphere of radiation

Oppenheimer and Volkoff write $e^{-\lambda} = 1 - 2u/r$, their (8), and for the Fermi gas at $t_0 \to \infty$ find the exact solution $e^t = 3/7r^2$, $u = 3r/14$, their (22).
With their (20), $du/dr = \tfrac12 r^2 e^t = 4\pi\rho r^2$, that is $8\pi\rho = 3/7r^2$, and $e^\lambda = 7/4$.
Their footnote 15 says it is a limiting form of Tolman's solutions V and VI.

Tolman's solution VI at $n = 1/2$, his (8.1), is $e^\lambda = 7/4$ and $e^\nu = (Ar^{1/2} - Br^{3/2})^2$, which at $B = 0$ is $e^\nu = A^2r$.
His solution V at $n = 1/2$, his (7.1), is $e^\lambda = 7/(4 - 7(r/R)^{7/3})$ and $e^\nu = B^2r$, which as $R \to \infty$ is the same metric.
He names it at the end of section 8: the limiting form "might be called the blackbody radiation solution".

The chart `areal` writes $e^\nu = r/a$ with $a$ a length, so that $t$ carries a time: $ds^2 = -(r/a)c^2dt^2 + \tfrac74dr^2 + r^2d\Omega^2$.
`misner_zapolsky_check` in `print_charts.py` holds it to $8\pi\rho = 3/7r^2$ and $p = \rho/3$, to $g^{rr} = 1 - 2(3r/14)/r$, to a vanishing Ricci scalar, and to the homothety $\xi = r\partial_r + \tfrac12 ct\,\partial_{ct}$, $\mathcal{L}_\xi g = 2g$.

## Step 3. The Tolman V star

Tolman's (7.1) to (7.9): the pressure $8\pi p = 1/7r^2 - (2/R^2)(r/R)^{1/3}$ vanishes at $r_b = R/14^{3/7}$, his (7.6), the vacuum outside fixes $B^2 = 14^{3/7}/2R$, his (7.8), and the mass is $m = r_b/4$, his (7.9).
In the radius of the star, $(r/R)^{7/3} = (r/r_b)^{7/3}/14$ and $B^2 = 1/2r_b$, so
$$g_{tt} = -\frac{r}{2r_b}, \qquad g^{rr} = Z = \frac47 - \frac{1}{14}\left(\frac{r}{r_b}\right)^{7/3}.$$
The chart `tolman_v` takes $r_b$ for its one parameter and names $Z$.

$Z$ holds a power $7/3$ of $r$, which the checker's zero test does not join to the other powers of $r$, so $Z$ is held as a function of $r$ with the declared slope $\partial_rZ = (7Z - 4)/(3r)$, `HELD` and `RATES` in `verify_metrics.py`, and every value of the chart is a rational function of $r$, $r_b$ and $Z$.
`misner_zapolsky_check` holds the slope to the definition, the density and pressure to his (7.2) and (7.3), the pressure to vanishing at $r_b$, and $g_{tt}$ and $g_{rr}$ there to Schwarzschild's published values at $r_s = r_b/2$.

## Step 4. The Tolman VI star

Tolman's (8.1) to (8.9): $8\pi\rho = 3/7r^2$, the pressure $8\pi p = (1/7r^2)(1 - 9(B/A)r)/(1 - (B/A)r)$ vanishes at $r_b = A/9B$, his (8.6), the vacuum outside fixes $A^2 = (3^6/2^4\cdot7)(B/A)$, his (8.8), and the mass is $m = 3r_b/14$, his (8.9).
In the radius of the star, $B/A = 1/9r_b$ and $A^2 = 81/112r_b$, so
$$g_{tt} = -\frac{r}{112r_b}\left(9 - \frac{r}{r_b}\right)^2, \qquad g_{rr} = \frac74.$$
`misner_zapolsky_check` holds the chart `tolman_vi` to his (8.2) and (8.3), the pressure to vanishing at $r_b$, and $g_{tt}$ and $g_{rr}$ there to Schwarzschild's published values at $r_s = 3r_b/7$.
`_tools/test_misner_zapolsky.py` holds his equation of state (8.5) as well.

## Step 5. Tolman's exponent

Tolman's solution V, his (4.5), with $n$ free and $R \to \infty$: $e^\nu = B^2r^{2n}$, $e^\lambda = 1 + 2n - n^2$, $8\pi\rho = (2n - n^2)/((1 + 2n - n^2)r^2)$ and $8\pi p = n^2/((1 + 2n - n^2)r^2)$, so $p/\rho = n/(2 - n)$.
Chavanis (2008) writes the same family in $q = p/\epsilon$, his (27) and (28): $e^\nu = Ar^{4q/(1+q)}$ and $e^{-\lambda} = 1 - 4q/((1+q)^2 + 4q)$, which is Tolman's at $n = 2q/(1 + q)$.
The chart `power_law` writes $e^\nu = (r/a)^{2n}$.
`misner_zapolsky_check` holds it to Tolman's density and pressure and to being the areal chart at $n = 1/2$, and the printer writes each value as a rational function of $n$ and $r$ times a power $2n$ of $r/a$, through `named_powers`.

## Step 6. How a star of finite central density settles onto the sphere

With $u = m/r$ and $v = 4\pi\rho r^2$ against $s = \ln r$, hydrostatic equilibrium for $p = \rho/3$ is
$$\frac{du}{ds} = v - u, \qquad \frac{dv}{ds} = 2v - \frac{4v\left(u + v/3\right)}{1 - 2u}.$$
The sphere of radiation is the fixed point $u = v = 3/14$.
The matrix of the linearised flow there has trace $-3/2$ and determinant $7/2$, so a small departure goes as $e^{-3s/4}\cos(\tfrac{\sqrt{47}}{4}s + \delta)$, which is Chavanis's (2002) expansion (134) at $q = 1/3$.
Integrated from a regular centre, $2m/r$ rises to $0.4926$ at $\alpha = \sqrt{16\pi\rho_c}\,r = 4.697$, where the central density is $22.40$ times the local one, and then oscillates about $3/7$: Chavanis's (2008) $\chi_c = 0.493$, $\alpha_c = 4.7$ and density contrast $22.4$.
`TheApproach` in `_tools/test_misner_zapolsky.py` holds all of it in plain Python.

## Step 7. What the diagrams draw

The spacetime diagrams draw the plane of $t$ and $r$ of each chart: the sphere of radiation at $a = 1$, each star at $r_b = 1$, and Tolman's exponent at $n = 1$, $a = 1$.
A radial ray has $dr/d(ct) = \sqrt{-g_{tt}/g_{rr}}$, which for the sphere of radiation is $\sqrt{4r/7a}$, so a ray from $r$ reaches the centre at $ct = \sqrt{7ar}$.
For Tolman's exponent the time to the centre is $\sqrt{1 + 2n - n^2}\,a^nr^{1-n}/(1 - n)$, infinite at $n = 1$, while the affine parameter $\int\sqrt{-g_{tt}g_{rr}}\,dr$ goes as $r^{n+1}$ and is finite.

The conformal diagrams use the tortoise coordinate $r_* = \int\sqrt{g_{rr}/(-g_{tt})}\,dr$ and $p, q = \arctan((ct \mp r_*)/\ell)$.
For the sphere of radiation $r_* = \sqrt{7ar}$, zero at the centre, so the diagram is Minkowski's triangle with a timelike singular centre.
Each star is joined to Schwarzschild's published exterior, with $g_{tt}$ and $g_{rr}$ checked continuous at $r_b$; inside solution VI $r_* = \tfrac{14}{3}r_b\ln((3 + \sqrt{r/r_b})/(3 - \sqrt{r/r_b}))$, checked against the quadrature.
At $n = 1$, $r_* = \sqrt2\,a\ln(r/a)$ covers the whole line, the diagram is Minkowski's diamond, and the centre is on its two null edges on the left, checked to be reached by an ingoing ray at the affine parameter $a/\sqrt2$ from $r = a$.

The embedding diagrams draw the equator at one moment of $t$.
With $g_{rr}$ constant the surface is a cone, $dz/dr = \sqrt{g_{rr} - 1}$: $\sqrt3/2$ for radiation, a flat sheet with a wedge of $2\pi(1 - 2/\sqrt7)$ cut out, and $1$ for the stiffest fluid.
The Tolman V star climbs at $\sqrt{(1 - Z)/Z}$, from $\sqrt3/2$ at the apex to $1$ at its surface, and the Tolman VI star is the cone of radiation cut off at $r_b$.
Each meets Flamm's paraboloid of its own $r_s$ with one tangent, since $g_{rr}$ is continuous there, which `embedding.misner_zapolsky` checks.
