# The three-brane and its throat

The extreme three-brane of type IIB supergravity in ten dimensions, and anti-de Sitter space of five dimensions times a 5-sphere, its throat.
The five charts are written by `_tools/derivations/print_charts.py --metric three_brane_throat`, which took 16 seconds for the slowest on 2 October 2026, and `verify_metrics.py --system three_brane_throat/<chart>` checks each in a few seconds.

## Step 1. The papers and the source of each chart

Every paper of the candidate entry was verified on Crossref and INSPIRE before anything was built.
Horowitz and Strominger (1991) and Duff and Lu (1991) are not on arXiv and were read in full from the KEK scans of their preprints, UCSBTH-91-06 (`lib-extopc.kek.jp/preprints/PDF/2000/0034/0034882.pdf`) and CTP/TAMU-29/91 (`.../1991/9107/9107505.pdf`).
Gibbons and Townsend (1993), Gibbons, Horowitz, and Townsend (1995), and Maldacena (1998) were read from arXiv, hep-th/9307049, hep-th/9410073, and hep-th/9711200.
Duff and Stelle (1991), Güven (1992), Dabholkar, Gibbons, Harvey, and Ruiz Ruiz (1990), and Callan, Harvey, and Strominger (1991) were verified by their records and abstracts alone, and the History says of them only what Duff and Lu's introduction and Gibbons and Townsend's letter say.

- `isotropic`: Duff and Lu's (2.2) with (2.20) and (3.4), $ds^2 = e^{2A}\eta_{\mu\nu}dx^\mu dx^\nu + e^{-2A}\delta_{mn}dy^m dy^n$ with $e^{-4A} = 1 + |Q|/r^4$; their $|Q|$ is written $L^4$, their $r$ is written $\rho$, and $e^{-4A}$ is written $H$.
- `areal`: Horowitz and Strominger's (36) at $r_+ = r_- = L$, which is Gibbons and Townsend's (3) with $\gamma_x = 1/2$ and Gibbons, Horowitz, and Townsend's (3.12) at $d = 7$, $p = 3$. Their (4.1), $r^4 = \rho^4 + L^4$, carries it to the isotropic chart.
- `horizon`: Gibbons, Horowitz, and Townsend's coordinate of (2.14), the power $1/(p + 1)$ of $1 - (\mu/r)^{d-3}$, which they state for $d = 4$ and carry to any $d$ below (3.12); at $d = 7$, $p = 3$ it is $w = (1 - L^4/r^4)^{1/4}$, without their constant factor. Every component is even in $w$, their reflection (2.22).
- `throat`: Maldacena's (2.3), with his $U = r/\alpha'$ written back as $r$ and $L^4 = 4\pi gN\alpha'^2$, which is his (7.3).
- `throat_proper`: Gibbons and Townsend's (5), $e^{2\rho/a}(-dt'^2 + dx' \cdot dx') + d\rho^2 + a^2d\Omega^2$, with their $\rho$ written $\sigma$ and the sign of $\sigma$ chosen so that the horizon is at $\sigma \to -\infty$.

Duff and Lu's chart with the harmonic function left free in the six Cartesian coordinates across the brane, the chart of many parallel three-branes, is not printed: its curvature holds every second derivative of $H$ in six variables.
No chart that is regular across the horizon in every component is in the papers read; Gibbons, Horowitz, and Townsend's argument is that $r$ is analytic in $w$ while $t$, $x$, $y$, and $z$ end there, as the Poincaré coordinates of anti-de Sitter space do.

## Step 2. The root of $1 - w^4$

The checker's canonical form writes a root whose radicand leads with a negative term on its other branch, as `wahlquist.md` Step 5 and `tolman_vii.md` Step 4 record.
So the horizon chart holds the areal radius $r = L(1 - w^4)^{-1/4}$ as a function of $w$, `HELD`, with the slope $\partial_wr = w^3r/(1 - w^4)$ declared in `RATES`, and its line element is written in $w$ and $r$ alone, $w^2(-c^2dt^2 + \ldots) + r^2(dw^2/(w^2(1 - w^4)^2) + d\Omega_5^2)$.
$L$ stands in neither: it is the constant of integration of the slope, so $w$ and $r$ have no relation that a value could hold, and every value is a rational function of the two.

## Step 3. What the script checks

`three_brane_check` holds each chart to Duff and Lu's (3.1), $R_{MN} = F_{MPQRS}F_N{}^{PQRS}/96$, for the 5-form of the potential $A_{txyz}$ equal to the square of the warp factor, their (2.3) with (2.20): the form's part on the sphere is its dual, is closed, and has the flux $4L^4$ times the volume of the unit 5-sphere; the Ricci scalar vanishes, and the mixed Ricci tensor is $-k$ on the first five dimensions and $+k$ on the sphere, with $k = 4L^8/(\rho^4 + L^4)^{5/2}$.
The areal and horizon charts are the isotropic chart pulled back, the horizon chart is even in $w$, and the inversion $\rho \to L^2/\rho$ multiplies the isotropic metric by $L^2/\rho^2$, Gibbons and Townsend's (7).
The throat is the limit $\lambda \to 0$ of the isotropic chart at $\rho = \lambda r$ with $t$, $x$, $y$, and $z$ divided by $\lambda$, has no Weyl tensor, and has $R^M{}_N = \mp 4/L^2$; the proper distance chart is the throat pulled back along $r = Le^{\sigma/L}$.
A difference that sympy's `simplify` leaves standing is settled at three points to thirty digits.

## Step 4. The drawings

Everything is drawn at $L = 1$.
The spacetime diagrams are the plane of the time and the radial coordinate of each chart, and `null_rays.py --verify` holds every ray to $ct \mp \rho_*$, with $d\rho_*/d\rho = \sqrt{1 + L^4/\rho^4}$, which is $\rho\,{}_2F_1(-1/2, -1/4; 3/4; -L^4/\rho^4)$ less the constant $2\Gamma(3/4)^2/\sqrt{\pi}$ that makes it odd in $\rho$.
The conformal diagram is Gibbons, Horowitz, and Townsend's Figure 1, a tower of exteriors alternately on the right and on the left: the region behind a horizon is the exterior under $w \to -w$, where $\rho_*$ continues as an odd function, so the advanced time of the exterior is the retarded time of the mirror region.
The throat alone is drawn as the Poincaré wedge of the strip of anti-de Sitter space of two dimensions, as Bertotti and Robinson's is.
The embedding diagram has the surface of $\rho$ and $\phi$ at one moment, with circumference radius $(\rho^4 + L^4)^{1/4}$ and $dz/d\rho = \sqrt{2\rho^4 + L^4}\,L/(\rho(\rho^4 + L^4)^{3/4})$, a plane that opens into an infinitely long throat, and the throat alone, the cylinder $d\sigma^2 + L^2d\phi^2$.
The spacetime is static, so there is no movie and no stack, and no new kind of diagram.
