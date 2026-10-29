# The regular black hole from pure gravity

Bueno, Cano and Hennigar, "Regular black holes from pure gravity", Phys. Lett. B 861 (2025) 139260, arXiv:2403.04827.
`print_charts.py` writes all three charts; `regular_black_hole_pure_gravity_check` holds them to what follows.

## Step 1: the field equation

For the theory $I = \frac{1}{16\pi G}\int d^Dx\sqrt{|g|}\,[R + \sum_{n \ge 2}\alpha_n\mathcal{Z}_n]$ of quasi-topological densities $\mathcal{Z}_n$, and the ansatz $ds^2 = -N^2f\,dt^2 + dr^2/f + r^2d\Omega_{D-2}^2$, the field equations are $N' = 0$ and $[r^{D-1}h(\psi)]' = 0$ with $h(\psi) = \psi + \sum_n\alpha_n\psi^n$ and $\psi = (1 - f)/r^2$, their (3) to (5).
So $h(\psi) = m/r^{D-1}$ with $m$ fixed by the mass, $M = (D-2)\Omega_{D-2}m/(16\pi G)$, their (8) with $h(\psi_+)r_+^{D-1} = m$.

## Step 2: the geometric series

At $\alpha_n = \alpha^{n-1}$, $h = \psi/(1 - \alpha\psi)$ and $f = 1 - mr^2/(r^{D-1} + \alpha m)$, their (7) and Table 1.
At $D = 5$, $\Omega_3 = 2\pi^2$, so $m = 8GM/(3\pi c^2)$, an area, which the page writes $r_0^2$, and $\alpha$, an area too, the page writes $\ell^2$:
$f = 1 - r_0^2r^2/(r^4 + \ell^2r_0^2)$.
At $D = 4$ the same function is Hayward's $F = 1 - 2m_Hr^2/(r^3 + 2m_H\ell^2)$ with $m = 2m_H$, $\alpha = \ell^2$; the check reads Hayward's published $g^{rr}$ and compares.

## Step 3: horizons, centre, far away

$r^4f = r^4 - r_0^2r^2 + \ell^2r_0^2$, so $r_\pm^2 = \tfrac{1}{2}(r_0^2 \pm r_0\sqrt{r_0^2 - 4\ell^2})$: two horizons for $r_0 > 2\ell$, one at $r^2 = 2\ell^2$ for $r_0 = 2\ell$, none below.
The surface gravity is $f'(r_+)/2 = (r_+^2 - 2\ell^2)/r_+^3$, which is $2\pi$ times their $T = (r_+^2 - 2\alpha)/(2\pi r_+^3)$ at $D = 5$, their (10).
Near the centre $f = 1 - r^2/\ell^2 + O(r^6)$: de Sitter space of radius $\ell$, where the printed $K = 40/\ell^4$ and $R = 20/\ell^2$.
Far away $f = 1 - r_0^2/r^2 + O(r^{-6})$, Tangherlini's with $r_h = r_0$, and at $\ell = 0$ the static chart is Tangherlini's published chart.

## Step 4: the charts

The static chart is the ansatz itself at $N = 1$, as Tangherlini's of 1963 is.
The Eddington-Finkelstein charts are Eddington's (1924) and Finkelstein's (1958) construction on the tortoise coordinate $dr_*/dr = 1/f$, with $v = ct + r_*$ and $u = ct - r_*$; the check pulls the static chart back along each.

## Step 5: the drawings

All drawings take $r_0 = 1$ and $\ell = 12/25$, where $r^4 - r^2 + \ell^2 = (r^2 - 16/25)(r^2 - 9/25)$, so $r_+ = 4/5$ and $r_- = 3/5$.
Then $1/f = 1 + r^2/((r^2 - r_+^2)(r^2 - r_-^2))$ and $r_* = r + \tfrac{10}{7}\ln|(r - \tfrac45)/(r + \tfrac45)| - \tfrac{15}{14}\ln|(r - \tfrac35)/(r + \tfrac35)|$, with $\kappa_+ = 7/20$ and $\kappa_- = 7/15$.
The embedding slope $dz/dr = r/\sqrt{(r^2 - r_+^2)(r^2 - r_-^2)}$ integrates to $\tfrac12\operatorname{arcosh}((2r^2 - r_+^2 - r_-^2)/(r_+^2 - r_-^2))$ outside $r_+$ and its mirror inside $r_-$.
$R$ vanishes at $r^4 = 5\ell^2r_0^2/3$, between the horizons, and is least, $-1/(4\ell^2)$, at $r^4 = 3\ell^2r_0^2$, outside $r_+$.
