# Bonnor's stars of charged dust

Charged dust at rest whose charge density equals its mass density, $\sigma = \pm\sqrt{4\pi\epsilon_0G}\,\rho$, in Majumdar and Papapetrou's metric $-U^{-2}c^2dt^2 + U^2$ times flat space.
Inside the dust the flat Laplacian of $U$ is $-4\pi G\rho\,U^3/c^2$, and outside it $U$ is harmonic.
`print_charts.bonnor_charged_dust` writes the eight charts, `bonnor_charged_dust_check` holds each to the Einstein-Maxwell equations with the dust for a source, and `_tools/test_bonnor_charged_dust.py` holds the published files to the same facts.

## Step 1. The sources, each read before anything was built

- A. Das, "A class of exact solutions of certain classical field equations in general relativity", Proc. R. Soc. A 267, 1 (1962), doi:10.1098/rspa.1962.0079. Crossref's record and abstract. His first name, Anadi, is from the Mathematics Genealogy Project, id 17959.
- W. B. Bonnor, "The Equilibrium of a Charged Sphere", MNRAS 129, 443 (1965), doi:10.1093/mnras/129.6.443. Read in full from the ADS scan.
- W. B. Bonnor and S. B. P. Wickramasuriya, "Are Very Large Gravitational Redshifts Possible?", MNRAS 170, 643 (1975), doi:10.1093/mnras/170.3.643. Read in full from the ADS scan.
- J. P. S. Lemos and E. J. Weinberg, "Quasiblack holes from extremal charged dust", Phys. Rev. D 69, 104004 (2004), arXiv:gr-qc/0311051. Read in full.
- J. P. S. Lemos and V. T. Zanchin, "Bonnor stars in d spacetime dimensions", Phys. Rev. D 77, 064003 (2008), arXiv:0802.0530. Abstract.

## Step 2. The charts and where each comes from

| chart | source | potential |
| --- | --- | --- |
| `harmonic` | Bonnor and Wickramasuriya's (2.10), Lemos and Weinberg's (2.8) | $U(x, y, z)$ left free |
| `sphere_1965` | Bonnor's (2.5) and (3.10), his $f = 1/U$ and his $a$ written $r_0$ | $U = (r_0 + m)^{3/2}/\sqrt{r_0^3 + mr^2}$ |
| `sphere_1975` | Bonnor and Wickramasuriya's (3.2) | $U = 1 + (m/2r_0)(3 - r^2/r_0^2)$ |
| `exterior` | Bonnor's (3.2), Bonnor and Wickramasuriya's (3.1) | $U = 1 + m/r$ |
| `exterior_areal` | Bonnor's (3.5) at $e^2 = m^2$, Bonnor and Wickramasuriya's (3.8), their $\bar r$ written $R$ | $R = r + m$ |
| `spheroid_interior` | Bonnor and Wickramasuriya's (4.1) and (4.3) | $U = 1 + (m/a)(\arctan(1/\sinh u_0) + (u_0^4 - u^4)/(4u_0^3\cosh u_0))$ |
| `spheroid_exterior` | Bonnor and Wickramasuriya's (4.1) and (4.2) | $U = 1 + (m/a)\arctan(1/\sinh u)$ |
| `quasi_black_hole` | Lemos and Weinberg's (3.1), their $R$, $q$ and $c$ written $r$, $m$ and $b$ | $U = 1 + m/\sqrt{r^2 + b^2}$ |

Lemos and Weinberg's $c$ is written $b$ because $c$ is the speed of light on every page.
Bonnor and Wickramasuriya's interior in the areal radius, their (3.9), holds the isotropic radius as a root of a cubic, $\bar r = rU(r)$, so it is not printed; the areal chart is given outside the star, where $R = r + m$.

## Step 3. A named potential is held with its slope

Each chart that names $U$ lists it in `HELD` and its slope in `RATES` of `verify_metrics.py`: $-mrU^3/(r_0 + m)^3$, $-mr/r_0^3$, $-mu^3/(au_0^3\cosh u_0)$, $-m/(a\cosh u)$ and $-r(U - 1)^3/m^2$.
The reader checks each slope against the definition, and every tensor is then a rational function of $U$ and the coordinates, as the line element is.
The sphere of 1965 holds $r_0$ only through $(r_0 + m)^3$, which is printed as that sum, and the spheroid's charts are printed in $\sinh u$, $\cosh u$ and $\cosh u_0$: the exponential of $u_0$ is written $\cosh u_0 + \sinh u_0$ and the $\sinh u_0$ cancels, which the printer checks.

## Step 4. What the check holds

In units $G = c = 4\pi\epsilon_0 = 1$, with $A = U^{-1}dt$, $u_a = -U^{-1}dt$ and $4\pi\rho = -\nabla^2U/U^3$, the Einstein tensor is $8\pi\rho\,u_au_b$ plus the Maxwell stress, and Maxwell's equations have the current $\rho u^b$, in every chart.
The densities are Bonnor's (3.11), Bonnor and Wickramasuriya's (3.3) and (I.1) and Lemos and Weinberg's (3.2), each positive, the first falling outward and the second rising.
$U$ and its slope are continuous across each surface.
The sphere of 1975 has the central redshift $3m/2r_0$, and as $r_0 \to 0$ its proper radius tends to $4m/3$, its area to $4\pi m^2$ and its central density to $2c^2/(9\pi Gm^2)$, the three numbers on Bonnor and Wickramasuriya's page 646.
The spheroid's charts are the harmonic chart pulled back along their (4.10), its redshift on the disc is their (4.8), and $U \to 1 + m/r$ far away.
The isotropic exterior is the published single hole of `majumdar_papapetrou`, and the areal exterior is the isotropic one at $r = R - m$ and the published `rn_metric` at $r_s = 2m$, $r_q = m$.

## Step 5. The drawings

All are at $m = 1$: both spheres at $r_0 = 2m$, the spheroid at $a = m$ and $u_0 = 1$, the cloud at $b = m/2$, which is also the harmonic chart's declared $U$.
The spacetime diagrams are the plane of $t$ and the radius of each chart, the line through the centre of the sphere of 1975, and the spheroid's axis of symmetry; each ray keeps $ct \mp \int U^2$ times the flat length, in closed form but on the spheroid's axis, where it is a quadrature.
The conformal diagrams are Minkowski's triangle for each, the star a timelike tube.
The embedding diagrams are the equator of each sphere and of the cloud, $\rho = rU$ and $dz/dr = \sqrt{-r\,\partial_rU(2U + r\,\partial_rU)}$, and the spheroid's equatorial plane, whose disc inside the focal ring is flat.
