# The perturbed Friedmann universe in Bardeen's potentials

Why each chart is the one published, what the checker counts as an order when the small quantities are potentials of a cosmology, and what `print_charts.py`, `null_rays.py`, `embedding.py` and `curvature.py` check before they write.

## Step 1: the sources

James M. Bardeen, "Gauge-invariant cosmological perturbations", Physical Review D 22, 1882-1905 (1980), doi:10.1103/PhysRevD.22.1882, is the origin of the gauge-invariant variables, and it was not read: the publisher's copy is behind a paywall and no scan or reprint could be read, as `blocked-sources.md` of the task records.
Nothing here is taken from it.
Every equation is taken from three sources that were read whole at the level of their equations:

- V. F. Mukhanov, H. A. Feldman and R. H. Brandenberger, "Theory of cosmological perturbations", Physics Reports 215, 203-333 (1992), doi:10.1016/0370-1573(92)90044-Z, read in the published scan;
- Chung-Pei Ma and Edmund Bertschinger, "Cosmological Perturbation Theory in the Synchronous and Conformal Newtonian Gauges", The Astrophysical Journal 455, 7-25 (1995), doi:10.1086/176550, arXiv:astro-ph/9506072, read in its source on arXiv;
- James M. Bardeen, Paul J. Steinhardt and Michael S. Turner, "Spontaneous creation of almost scale-free density perturbations in an inflationary universe", Physical Review D 28, 679-693 (1983), doi:10.1103/PhysRevD.28.679, read in its preprint UPR-0202T, the KEK scan 0036508.

The History also rests on R. K. Sachs and A. M. Wolfe, The Astrophysical Journal 147, 73-90 (1967), doi:10.1086/148982, read in the scan of the NASA Astrophysics Data System, and on the abstract of the Planck Collaboration's cosmological parameters of 2018, arXiv:1807.06209.
E. Lifshitz's paper of 1946 could not be read, and is cited only as the origin of the synchronous gauge, which Ma and Bertschinger's section 1 and Mukhanov, Feldman and Brandenberger's chapter 1 both state.
Hideo Kodama and Misao Sasaki's review of 1984 could not be read, and nothing rests on it.

## Step 2: the charts

Mukhanov, Feldman and Brandenberger write the scalar perturbations of a Friedmann universe in the signature $(+,-,-,-)$ as their (2.10), $ds^2 = a^2(\eta)\{(1 + 2\phi)d\eta^2 - 2B_{|i}dx^id\eta - [(1 - 2\psi)\gamma_{ij} + 2E_{|ij}]dx^idx^j\}$, on the background (2.1) with $\gamma_{ij} = \delta_{ij}[1 + \tfrac{1}{4}\mathcal{K}(x^2 + y^2 + z^2)]^{-2}$, their (2.2).
Their gauge-invariant potentials are (3.13), $\Phi = \phi + a^{-1}[(B - E')a]'$ and $\Psi = \psi - (a'/a)(B - E')$, and below it they write "In his notation, $\Phi = \Phi_A$ and $\Psi = -\Phi_H$", his being Bardeen's of 1980.
Their longitudinal gauge is $B = E = 0$, which fixes the coordinates entirely, and its metric is (3.21), $ds^2 = a^2[(1 + 2\Phi)d\eta^2 - (1 - 2\Psi)\gamma_{ij}dx^idx^j]$.

Ma and Bertschinger write the same gauge, which they call the conformal Newtonian gauge, in their section 2 in the signature $(-,+,+,+)$ and with $c = 1$, $ds^2 = a^2(\tau)\{-(1 + 2\psi)d\tau^2 + (1 - 2\phi)dx^idx_i\}$ on flat space, and state there $\psi = \Phi_A$ and $\phi = -\Phi_H$.
Their $\psi$ is the $\Phi$ here and their $\phi$ the $\Psi$ here, so the two papers give Bardeen's potentials the same signs.
The page takes Mukhanov, Feldman and Brandenberger's letters, reversed into the signature $(-,+,+,+)$, which reverses the line element and changes no function:

- `longitudinal`: $ds^2 = a^2\left(-(1 + 2\Phi)d\eta^2 + (1 - 2\Psi)(dx^2 + dy^2 + dz^2)\right)$;
- `longitudinal_curved`: the same with the space of (2.2), $(1 - 2\Psi)(dx^2 + dy^2 + dz^2)/(1 + k(x^2 + y^2 + z^2)/4)^2$, the $\mathcal{K}$ of (2.2) written $k$ as the published charts of `frw` write their curvature;
- `synchronous`: $\phi = B = 0$ in (2.10), $ds^2 = a^2\left(-d\eta^2 + ((1 - 2\psi)\delta_{ij} + 2\partial_i\partial_jE)dx^idx^j\right)$, which is the synchronous gauge of Ma and Bertschinger's section 2 with $h_{ij} = -2\psi\,\delta_{ij} + 2\partial_i\partial_jE$, its scalar part.

The conformal time $\eta$ is a length, with $a\,d\eta = c\,dt$, so the components are taken in $\eta$ itself and no factor of $c$ appears; $x$, $y$ and $z$ are comoving lengths, $a$, $\Phi$, $\Psi$ and $\psi$ pure numbers, $E$ an area and $k$ an inverse area, as `DIMENSIONS` in `verify_metrics.py` declares.
The functions are left free, so no component assumes a field equation.

## Step 3: what the checker counts as an order

Each chart is kept to first order in its potentials: `ORDERS` declares `({"Phi": 1, "Psi": 1}, 1)` for the longitudinal charts and `({"psi": 1, "E": 1}, 1)` for the synchronous one, with no time named.
Until this spacetime a declared function could be counted as small only in a post-Newtonian system, which names a slow time along which each derivative counts one order more; `Reader` now counts it in any system, and where no time is named every derivative of a potential counts as the potential does, since a cosmological perturbation varies as fast in time as in space.
`Reader._scaled` is the whole of the change: with no slow time the power of a derivative along time is zero.
A system that names no time and declares no function counts constants alone, as before, and the post-Newtonian systems are unchanged; `_tools/test_bardeen_perturbations.py` holds both, and holds a function of all four coordinates added at second order to every component to move no tensor handed over.

## Step 4: what each chart is held to

`bardeen_perturbations_check` in `print_charts.py` runs before anything is written.

The two longitudinal charts are held to Mukhanov, Feldman and Brandenberger's Einstein tensor, the background (4.2) plus the perturbation (4.8) at $B = E = 0$, written on the space of (2.2) with its Laplacian and its covariant Hessian.
Reversing the signature keeps $G_{\mu\nu}$, since the Ricci tensor is unchanged and the metric and the Ricci scalar both change sign, and so reverses $G^\mu{}_\nu$: every mixed component of each chart is minus theirs, to first order.
$G^i{}_0$ is not in (4.8); at first order $G^i{}_0 = -\gamma^{ij}G^0{}_j$, since both are the one first order $G_{0i}$ raised by the background metric.
In flat space the four equations are Ma and Bertschinger's Einstein equations of the conformal Newtonian gauge in their section 4, in Fourier space, with $\delta T^0{}_0 = -\delta\rho$ and $\nabla^2 \to -k^2$.

The synchronous chart is held to the longitudinal one pulled back along Mukhanov, Feldman and Brandenberger's (3.22), $\eta_l = \eta_s + \partial_\eta E$ and $x_l^i = x_s^i + \partial_iE$, with their (3.23), $\Phi = -\partial_\eta^2E - (a'/a)\partial_\eta E$ and $\Psi = \psi + (a'/a)\partial_\eta E$, the scale factor moved to $a + a'\partial_\eta E$ and the potentials, which move at second order only, left where they are.
The two metrics agree to first order in every component, so Bardeen's potentials of the synchronous chart are $\Phi_A = -\partial_\eta^2E - (a'/a)\partial_\eta E$ and $\Phi_H = -\psi - (a'/a)\partial_\eta E$, which its convention states.

## Step 5: what the components say

Every value is printed by `bardeen_perturbations_grouped` as a sum over the products of the scale factor and its derivatives, each times the sum of its terms in the potentials, and by `bardeen_perturbations_order` with the background first.
The Friedmann universe's own term then stands with its first order correction beside it: $R = 6(1 - 2\Phi)a''/a^3 - \dots$ and $K = 12(1 - 4\Phi)a''^2/a^6 - \dots$, and at $\Phi = \Psi = 0$ each value is the published conformal chart of `frw` at $k = 0$, or at its $k$ for the curved chart.
$G_{xy} = \partial_x\partial_y(\Psi - \Phi)$, so a stress with no off diagonal part gives $\Phi = \Psi$, and $G_{\eta\eta} = 3a'^2/a^2 - 6\partial_\eta\Psi\,a'/a + 2\nabla^2\Psi$ is the time equation, Poisson's with the expansion's two terms.

## Step 6: the drawings

Every drawing takes the width $L$ of one well as its unit of length, and a Gaussian well of depth $1/10$ held still in comoving coordinates, $\Phi = \Psi = -\tfrac{1}{10}e^{-(x^2 + y^2 + z^2)/L^2}$, `nr.BP_WELL`.
In the flat chart the universe is dust, $a = \eta^2/L^2$, Ma and Bertschinger's $a \propto \tau^2$ of the matter era, and a constant potential is Mukhanov, Feldman and Brandenberger's growing mode, their (17.1) and (17.2), $\Phi = A(\mathbf{x})$.
With $\Phi = \Psi$ constant and $a = \eta^2$ the published pressure perturbation vanishes, since $2\mathcal{H}' + \mathcal{H}^2 = 0$ for $\mathcal{H} = 2/\eta$, so the declared universe is a solution for dust, which `curvature.py` checks on the published $G^x{}_x$.
The curved chart is drawn at $k = 1/4L^2$, a closed universe of dust, $a = 1 - \cos(\eta/2L)$, with the same well held still; there a constant potential is no solution, and the view says only what the metric does with the well declared.
The synchronous chart is drawn as the flat universe and its well carried along Mukhanov, Feldman and Brandenberger's (3.25) and (3.26): $\int a\Phi\,d\eta = \Phi\eta^3/3$, so $E = -\int d\eta\,a^{-1}\int a\Phi\,d\eta = -\eta^2\Phi/6$ and $\psi = \Phi + (a'/a^2)\int a\Phi\,d\eta = \tfrac{5}{3}\Phi$; the convention's formulas then give $\Phi_A = \Phi/3 + 2\Phi/3 = \Phi$ and $\Phi_H = -\tfrac{5}{3}\Phi + \tfrac{2}{3}\Phi = -\Phi$.

The spacetime diagrams draw the plane of $\eta$ and $x$ through the centre of the well, $y = z = 0$, from $\eta = 0$ to $3L$ and $x$ from $-1.5L$ to $1.5L$.
The two longitudinal views are checked by `null_rays.py --verify` against $\eta \pm x_*(x)$, with $x_*$ the quadrature `_bp_well` of $\sqrt{(1 - 2\Psi)/(1 + 2\Phi)}/\Omega$; on 23 October 2026 both families held to $3 \times 10^{-12}$ and $1.4 \times 10^{-11}$.
The synchronous chart has no such form, since its metric on the plane depends on $\eta$ and $x$ together; its coordinates are first order alone, and $g_{xx}$ at the centre, $1 + 1/3 - \eta^2/15$ in units of $a^2$, stays above $0.7$ on the plane drawn.

The embedding diagram draws the plane $z = 0$ at $\eta = L$, where $a = 1$, as a surface of revolution about the centre of the well, $\rho = s\sqrt{1 - 2\Psi}$, out to $s = 3L$; it rises by $0.62\,L$, and the rise goes as the square root of the depth of the well.
Its moment is the line $\eta = L$ on the flat view, and on the synchronous view the curve $\eta(1 + e^{-x^2}/30) = 1$, the same moment carried along $\eta \to \eta + \partial_\eta E$; the curved view draws another spacetime, which `slices.HIDDEN` says.

The curvature figure draws $K$ and $R$ along the same line at the same moment, beside the published scalars of `frw` at $k = 0$ and $a = \eta^2$: $R = 12 + (18 - 4x^2)e^{-x^2}/5$ and $K = 240 + (144 - 32x^2)e^{-x^2}$ in units of $L$.
For dust $R = 8\pi G\rho/c^2$, so $R/12$ is the density over its mean, $1.3$ at the centre, which equals the published $G^\eta{}_\eta$ over Friedmann's.
