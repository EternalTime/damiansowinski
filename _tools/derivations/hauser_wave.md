# Hauser's twisting type N wave

## Sources and scope

The founding record is I. Hauser, Physical Review Letters 33 (1974) 1112-1113, DOI 10.1103/PhysRevLett.33.1112, received 26 August 1974.
Its body has not been read: the APS full text needs a subscription and no open copy exists.
Its erratum, Physical Review Letters 33 (1974) 1525, is read in full and corrects Eq. (15) to $f'' + 3f/[16(1 + y^2)] = 0$.
The readable metric source is H. Godazgar, M. Godazgar, R. Monteiro, D. Peinador Veiga and C. N. Pope, Physical Review Letters 126 (2021) 101103, arXiv:2010.02925, eq. (24) for the twisting algebraically special vacua and eq. (29) for Hauser's functions.
F. J. Chinea and F. Navarro-Lérida, Classical and Quantum Gravity 17 (2000) 4587, rederive the field and end on the same equation for f (their eq. 75).
Firstmate approved this source base on 2026-10-30.

## Step 1. The metric of Godazgar et al

In signature (-+++), with $\partial = \partial_\zeta - L\partial_u$,
$ds^2 = -2\Omega\left(dr + W d\zeta + \bar W d\bar\zeta + H\Omega\right) + 2(r^2 + \Sigma^2)P^{-2}d\zeta d\bar\zeta$, $\Omega = du + L d\zeta + \bar L d\bar\zeta$,
$2i\Sigma = P^2(\bar\partial L - \partial\bar L)$, $W = -(r + i\Sigma)\partial_u L + i\partial\Sigma$, $H = K/2 - r\partial_u\ln P$, $K = 2P^2\,\mathrm{Re}\,\partial(\bar\partial\ln P - \partial_u\bar L)$.
Hauser's field is $L = 2i(\zeta + \bar\zeta)$, $P = (\zeta + \bar\zeta)^{3/2}f(t)$, $t = u/(\zeta + \bar\zeta)^2$, with $16(1 + t^2)f'' + 3f = 0$.

## Step 2. Real coordinates

Put $\zeta = x + iy$, so $\zeta + \bar\zeta = 2x$, $L = 4ix$ and $t = u/(4x^2)$.
Then $\bar\partial L - \partial\bar L = 4i$, so $\Sigma = 2P^2 = 16x^3f^2$ exactly.
The twist of the rays, $\Sigma/(r^2 + \Sigma^2)$, is therefore $2P^2/(r^2 + 4P^4)$, nonzero wherever $f$ is.
With $f$ real, $W = 32x^2ff' + 8ix^2f(3f - 4tf')$, so $W d\zeta + \bar W d\bar\zeta = 64x^2ff'\,dx - 16x^2f(3f - 4tf')\,dy$.
$L d\zeta + \bar L d\bar\zeta = -8x\,dy$, and $du = 4x^2dt + 8xt\,dx$, so $\Omega = 4x^2dt + 8xt\,dx - 8x\,dy$.
$K = 2x\left(8(1 + t^2)ff'' - 8(1 + t^2)f'^2 + 12tff' - 3f^2\right)$, and with Hauser's equation used once, $8(1 + t^2)ff'' = -3f^2/2$, it is $x\left(24tff' - 16(1 + t^2)f'^2 - 9f^2\right)$.
$\partial_u\ln P = f'/(4x^2f)$, so $H = x\left(24tff' - 16(1 + t^2)f'^2 - 9f^2\right)/2 - rf'/(4x^2f)$.
The transverse factor is $2(r^2 + \Sigma^2)/P^2 = (r^2 + 256x^6f^4)/(4x^3f^2)$.
Every coordinate is a pure number, and the length $L$ multiplies the whole line element.
These were computed with sympy from Step 1 and are the printed chart.

## Step 3. Why t is a coordinate

The page's reader declares a function of coordinates, and $f$ depends on $u$ and $x$ only through $t$.
Taking $t$ itself as a coordinate makes $f$ a function of one coordinate, whose derivatives the reader writes as primes.

## Step 4. The vacuum

`hauser_wave_check` in `print_charts.py` builds the Ricci tensor with $f$ free, and at three points of rational coordinates sets $f$ and $f'$ to rational values and takes $f''$ to $f^{(5)}$ from Hauser's equation.
Every component is then exactly zero.
With the coefficient 3 changed to 4 they are not, so the check can fail.
The published tensors leave $f$ free, as Robinson and Trautman's leave $P$ and $H$ free, and `metric_tags.py` carries `vacuum` by `OVERRULED`.
An independent check in the original $(u, r, x, y)$ chart, with $f$ a 12-term series solution about the point, gave $|R_{ab}| < 10^{-13}$ against entries of order 20.

## Step 5. Type N and the symmetries

The same check finds $C_{abcr} = 0$ on Hauser's equation: $\partial_r$, the null tangent to the rays, is annihilated by the Weyl tensor, which is the condition for a fourfold principal null direction, and the Weyl tensor does not vanish.
$\partial_y$ is a Killing vector, since no component depends on $y$.
Under $x \to \lambda x$, $y \to \lambda y$, $r \to \lambda^3 r$ with $t$ fixed, $\Omega \to \lambda^2\Omega$, $dr + W d\zeta + \bar W d\bar\zeta \to \lambda^3(\ldots)$, $H \to \lambda H$ and the transverse term scales as $\lambda^5$, so $\xi = x\partial_x + y\partial_y + 3r\partial_r$ obeys $\mathcal{L}_\xi g = 5g$, which the check confirms at the same points.
That is a homothety, $\mathcal{L}_\xi g = 2\phi g$ with $\phi = 5/2$, the value Finley, Plebański and Przanowski and Chinea and Navarro-Lérida give for Hauser's field.

## Step 6. Hauser's function

Novoseller's even solution is $f_e = {}_2F_1(-1/8, -3/8; 1/2; -t^2)$, as Ivanov writes it, and the odd one is $f_o = t\,{}_2F_1(1/8, 3/8; 3/2; -t^2)$.
Both satisfy the equation to $10^{-15}$ at $t = 0.3, 2.5, -7$ (mpmath).
$f_e(0) = 1$, and its first zero is at $t_0 = 5.46728770730835$; where $f$ vanishes, $P$ does and the chart ends.
For large $|t|$ the equation's solutions grow as $|t|^{3/4}$ or $|t|^{1/4}$, the roots of $16\alpha(\alpha - 1) + 3 = 0$.
The drawings take $f = f_e$ on $|t| < t_0$.

## Step 7. No conformal diagram

Hauser's field is not asymptotically flat (Sommers and Walker 1976), it is singular on hypersurfaces (Chudecki and Przanowski 2018), and its one Killing vector is the translation $\partial_y$.
No plane of two coordinates carries its causal structure: the light rays of the plane of $t$ and $r$ leave it in $x$, and no compactification of the whole field is known.
So `conformal.py` names it in `NOT_DRAWN`.

## Step 8. No embedding diagram

No slice of this spacetime has an angle.
Its one symmetry is the translation $\partial_y$, and its wave surfaces of constant $t$ are null, since $\partial_r$ lies in them and is null.
The surface of $x$ and $y$ at fixed $t$ and $r$ is curved, but at $t = 0$, where $f = 1$ and $f' = 0$, its metric is $L^2\,\mathrm{diag}(T, T - 192x^3)$ with $T = (r^2 + 256x^6)/(4x^3)$, spacelike only for $r^2 > 768x^6$.
A surface of revolution drawn from it by reading $y$ as an angle would show a periodicity the spacetime does not have, a compactification presented as a fact of Hauser's wave.
So `embedding.py` names it in `NOT_DRAWN` and `build_mfs_data.py` in `WITHOUT_EMBEDDING`, decided with firstmate on 2026-10-30.

## Step 9. No curvature figure

The figure format stores the Kretschmann scalar and the Ricci scalar alone.
On Hauser's equation the Ricci tensor vanishes, so $R = 0$; and in a type N vacuum the Weyl tensor is $\Psi_4$ alone, which cannot pair with itself in the quadratic invariants, so $K = C_{abcd}C^{abcd} = 0$ as well.
Every invariant of the curvature and its first derivatives vanishes, and the first that does not, $9(2\rho)^4(2\bar\rho)^4\Psi_4^2\bar\Psi_4^2$, has second derivatives (Bičák and Pravda 1998, eq. 4.15), which the format cannot store.
The published $R$ and $K$ hold $f$ free and vanish exactly when Hauser's equation is used, which `_tools/test_hauser_wave.py` checks at points, so there is nothing to plot and the spacetime has no curvature file.

## Step 10. The spacetime diagram

The one view is the plane of $t$ and $r$ on the line $x = 1/2$, $y = 0$, where $4x^2 = 1$ and $t$ is the retarded coordinate $u$ itself, drawn at $L = 1$ with Novoseller's even $f$, the declared `hauser_f` of `null_rays.py`, whose derivatives sympy takes from Hauser's equation.
On that plane the metric is $-2\,dt\,dr - 2H\,dt^2$, so the rays of constant $t$ are null, and they are the geodesic rays with $r$ affine; the other null family obeys $dr/dt = -H$, and $\Gamma^x{}_{tt}$ turns it out of the plane, so the caption calls it null curves.
The window $|r| \le 2.5$, $|t + r| \le 2.5$ is square and keeps $|t| \le 5$, inside the first zero $t_0 = 5.467$ of $f$, where $P$ vanishes and the chart ends.
The two curves on which the published $g^{rr}$ vanishes are marked; the plane is carried onto itself by $(t, r) \to (-t, -r)$, since $f$ is even, $f'$ odd and so $H$ unchanged by it, and the two curves go into each other.
No closed form of the rays is known to check them against.
