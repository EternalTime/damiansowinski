# The Newman-Tamburino metrics

Ezra Newman and Louis Tamburino found in 1962 every vacuum metric whose principal null rays are geodesic, hypersurface orthogonal, sheared and expanding, J. Math. Phys. 3, 902 (1962).
Their metrics fall into a spherical class, $\rho^2 \neq \sigma\bar\sigma$, and a cylindrical class, $\rho^2 = \sigma\bar\sigma$.
The entry carries the spherical class and the limit of the cylindrical class, each as one chart, in the forms John Steele writes in section 2 of arXiv:gr-qc/0403031 (Class. Quantum Grav. 21, 3587 (2004)), with his signature $(+,-,-,-)$ flipped and his $u$ the chart's $x^0 = cu$.

Both charts are written by `_tools/derivations/print_charts.py --metric newman_tamburino`, and `verify_metrics.py --system newman_tamburino/spherical` and `--system newman_tamburino/cylindrical` check each.
The 1962 paper itself is behind the publisher's wall and was read only in its abstract; Steele's preprint and L. De Groote and N. Van den Bergh, arXiv:0802.1141 and arXiv:0804.1029, were read in full.
The literature review behind the entry, with every claim, query and quotation, is kept with the task's records, not in the repository.

## Step 1. The spherical class

Steele gives the spherical class as its inverse metric, in $u$, $r$ and $\zeta = x + iy$, with $A(u) = bu + c$, $a = A(\zeta\bar\zeta)^{1/2}$, $L = \tfrac12\ln\left((r + a)/(r - a)\right)$ and $R^2 = r^2 - a^2$:
its only components are $g^{12} = 1$, $g^{22}$, $g^{23}$, $g^{24}$, $g^{33} = -2s^3/(r + a)^2$ and $g^{44} = -2s^3/(r - a)^2$, with $s = \sqrt{x^2 + y^2}$.
With $g^{11} = g^{13} = g^{14} = 0$ the metric follows at once: $g_{12} = 1$, $g_{33} = 1/g^{33}$, $g_{44} = 1/g^{44}$, $g_{13} = -g^{23}/g^{33}$, $g_{14} = -g^{24}/g^{44}$ and $g_{11} = -g^{22} + (g^{23})^2/g^{33} + (g^{24})^2/g^{44}$.
The brackets of $g^{23}$ and $g^{24}$ collapse:
$g_{13} = (x/s^2)\left(L(r + a)^2/a - r\right)$, since $(r - 2a)(r + a) + 2a^2 = r(r - a)$, and $g_{14} = (y/s^2)\left(L(r - a)^2/a - r\right)$, since $(r + 2a)(r - a) + 2a^2 = r(r + a)$.
With $(r - a)^2s^2 + 4arx^2 = (r + a)^2x^2 + (r - a)^2y^2$,
$g_{11} = (2L/(a^2 s))\left(a\,r\,s^2 - L\left((r + a)^2x^2 + (r - a)^2y^2\right)\right)$.
Flipping the signature and writing $du$ as $c\,du$ gives the line element of the chart.
Steele's constant $c$ is a length, written $A_0$ here since $c$ is the speed of light, so $A = A_0 + b\,c\,u$ with $b$ a pure number.

`newman_tamburino_check` holds the chart's inverse metric to Steele's, slot by slot, before anything is written.
The covariant metric was first checked Ricci flat by hand at two random rational points to twenty digits (outside the repository), and the check asks again for every component of the Ricci tensor to vanish exactly.

## Step 2. The cylindrical limit

Steele writes the limit of the cylindrical class as
$ds^2 = 2\,du\,dr - x^{-2}\left(b + \ln(r^2x^4)\right)du^2 + 4(r/x)\,du\,dx - r^2dx^2 - x^2dy^2$,
noting that Stephani and his co-authors correct the $du^2$ coefficient of the 1962 paper in their (26.23).
His form is Ricci flat.
The limit De Groote and Van den Bergh print at the end of their section 3 of arXiv:0804.1029 is not, as printed: its $R_{ux}$ is $-2c_0^2\rho/x^3$, a slip in the preprint, so nothing is taken from it.

The logarithm needs a length, so $b = -\ln\ell^2$ and the $du^2$ coefficient is $x^{-2}\ln(r^2x^4/\ell^2)$.
$\ell$ is no parameter of the geometry: $u \to \mu u$, $r \to r/\mu$, $x \to \mu x$, $y \to y/\mu$ leaves every other term alone and sends $\ln(r^2x^4)$ to $\ln(\mu^2r^2x^4)$, so it moves $b$ by $\ln\mu^2$.
`newman_tamburino_check` holds the chart to Steele's line element with that $b$ and asks for a vanishing Ricci tensor.

## Step 3. Units

In the spherical class $r$ is an affine length and $x$ and $y$ are pure numbers, so $a$ and $A$ are lengths, $L$ and $s$ are pure, and $b$ is pure.
In the cylindrical limit $x$ is a pure number and $y$ a length, since $r^2dx^2$ and $x^2dy^2$ are areas, and $\ell$ is a length.
`DIMENSIONS` declares these.

## Step 4. The held names

Written out, $L$ is a logarithm of radicals, $\tfrac12\ln\left((r + A\sqrt{x^2 + y^2})/(r - A\sqrt{x^2 + y^2})\right)$, and $A = A_0 + b\,c\,u$ spreads every value into powers of $A_0$ and $b$.
The metric holds $u$ only through $a = A\,s$, so $a$ and $L$ are held as functions (`HELD`), with their slopes declared in `RATES`:
$\partial_u a = b\,c\,s$, $\partial_x a = a\,x/s^2$, $\partial_y a = a\,y/s^2$, $\partial_r L = -a/(r^2 - a^2)$, and along $u$, $x$ and $y$ the slope of $L$ is $r/(r^2 - a^2)$ times that of $a$.
The reader checks each slope against the definition, and every value is then a rational function of $a$, $L$, $r$, $x$, $y$, $s$ and $b$, with no $u$ left in it.
$L$ is transcendental over the rest, and $a$ stands in no relation with the coordinates once $u$ is gone, so a value that vanishes in the held names vanishes.
In the chart $x^0 = cu$ the slopes along $x^0$ carry no $c$, and `newman_tamburino` hands the printer the rates and the definition of $A$ with the $c$ taken off.
The spherical chart takes 222 seconds to print, most of it in factoring the scalars, and its check takes 36.
The Kretschmann scalar it prints, $K = -64r^2a^2\left(4(L + b)(x^4 - y^4) + 18x^2y^2 - 3(x^4 + y^4)\right)/\left((r^2 - a^2)^4s^2\right)$, agrees with one computed from the metric at a random rational point outside the repository to twenty digits.
It diverges where $r \to |a|$, the end of every ray.

## Step 5. Petrov type

At a random point of each chart the Weyl scalars in a null tetrad built on $\ell = \partial_r$ have $\Psi_0 = 0$, $\Psi_1 \neq 0$ and $I^3 \neq 27J^2$, with $I = \Psi_0\Psi_4 - 4\Psi_1\Psi_3 + 3\Psi_2^2$ and $J$ the determinant of the Weyl scalars' Hankel matrix: both are of type I, with $\partial_r$ a principal null direction that is not repeated.
In the cylindrical limit $\Psi_4$ vanishes in that tetrad as well, and its Kretschmann scalar is $K = 16\left(3 - \ln(r^2x^4/\ell^2)\right)/(r^4x^4)$, which vanishes on the surface $r\,x^2 = e^{3/2}\ell$ and changes sign there.
