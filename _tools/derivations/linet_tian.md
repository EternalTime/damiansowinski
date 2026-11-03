# Linet and Tian's cylinder

The static vacuum with cylindrical symmetry and a cosmological constant $\Lambda$, found by Bernard Linet, J. Math. Phys. 27, 1817 (1986), and Qingjun Tian, Phys. Rev. D 33, 3549 (1986).
Neither paper was read in full; every equation here is taken from the papers that were, and each chart names its source.

## The line element and its exponents

In the proper distance $\rho$ from the axis, da Silva, Wang, Paiva and Santos, Phys. Rev. D 61, 044003 (2000), their (1) and (2), and Žofka and Bičák, Class. Quantum Grav. 25, 015011 (2008), their (1.1) to (1.3), write

$$ds^2 = d\rho^2 + Q^{2/3}\left(-P^{2p_t}c^2dt^2 + P^{2p_z}dz^2 + P^{2p_\phi}d\phi^2/C^2\right),$$

with $Q = \sin(\sqrt{3\Lambda}\rho)/\sqrt{3\Lambda}$, $P = 2\tan(\sqrt{3\Lambda}\rho/2)/\sqrt{3\Lambda}$, $p_t = -(4\sigma^2 - 8\sigma + 1)/(3\Sigma)$, $p_z = (8\sigma^2 - 4\sigma - 1)/(3\Sigma)$, $p_\phi = -2(2\sigma^2 + 2\sigma - 1)/(3\Sigma)$ and $\Sigma = 4\sigma^2 - 2\sigma + 1$.
Since $Q/P = \cos^2(\sqrt{3\Lambda}\rho/2)$, $Q^{2/3}P^{2p} = \cos^{4/3}(\sqrt{3\Lambda}\rho/2)P^{2(p + 1/3)}$, and $p + 1/3$ is Levi-Civita's Kasner exponent along the same direction: $p_0 = 2\sigma/\Sigma$ along $t$, $p_2 = (1 - 2\sigma)/\Sigma$ along $\phi$ and $p_3 = 2\sigma(2\sigma - 1)/\Sigma$ along $z$, with $p_0 + p_2 + p_3 = p_0^2 + p_2^2 + p_3^2 = 1$.
That is the form of Griffiths and Podolský, Phys. Rev. D 81, 064015 (2010), their (3.6) and (3.7), and every chart of the page is written in it, with the exponents free on the page's parameters and held to Kasner's circle by `PARAMETER_RELATIONS` in `verify_metrics.py`, as Levi-Civita's own Kasner chart is.
The exponents $p_t$, $p_z$ and $p_\phi$ of the History add to zero because $p_0 + p_2 + p_3 = 1$.

The checker's reader takes a subscript only on a letter or after a command, so $p_\phi$ cannot be a parameter's name; the Kasner names $p_0$, $p_2$ and $p_3$, which the live Levi-Civita page already uses, carry the same content.

## The charts

- `toroidal`, $\Lambda > 0$: Griffiths and Podolský's (3.5) with their (3.2), whose $z$ is periodic and renamed $\psi$, with a conicity $B$ for it beside $C$; both are written below the line, as Žofka and Bičák write $C$, which is a constant rescaling of Griffiths and Podolský's $B$ and $C$.
- `cylindrical`, $\Lambda < 0$: the same form with $\tanh$ and $\sqrt{-3\Lambda}$, the substitution Žofka and Bičák state below their (1.3) and derive in their Appendix A; Tyc and Žofka, Class. Quantum Grav. 43, 105010 (2026), their (4).
- `gleiser`, either sign: Gleiser, Class. Quantum Grav. 38, 247001 (2021), his (1), from his arXiv:1810.07296, sec. III, whose $p_i$ are twice $p - 1/3$, so his $y^{1/3 + p_i/2}(1 - \Lambda y)^{1/3 - p_i/2}$ is $y^p(1 - \Lambda y)^{2/3 - p}$ in the Kasner exponent $p$.
  In $Y = y/(1 - \Lambda y)$ that is $(1 + \Lambda Y)^{-2/3}Y^p$, and his $dy^2/(3y(1 - \Lambda y))$ is $(1 + \Lambda Y)^2dy^2/(3Y)$; the page writes the chart so, with $Y$ held.
  With $\Lambda y = \sin^2(\sqrt{3\Lambda}\rho/2)$, $dy^2/(3y(1 - \Lambda y)) = d\rho^2$, and $P^2 = 4Y/3$, so each component is the proper distance chart's times a constant.
- `sarioglu_tekin`, $\Lambda < 0$: Sarıoğlu and Tekin, Phys. Rev. D 79, 087502 (2009), their (1) at $n = 3$, with $\ell = L$, $\Lambda = -3/L^2$ and their Kasner conditions (2); their (3), $r = r_0\cosh^{2/3}(3\rho/2L)$, takes it to the cylindrical chart, each component a constant times the other's.

Not built: Krasiński's own chart, J. Math. Phys. 16, 125 (1975), his (7.3) with $V = (x^2 - p_0)(x^2 - q_0)$ for real $p_0 < q_0$, his Type B, which Gleiser's (13) to (19) turn into the static form and his (23) writes out.
Its parameters are the roots $p_0$, $q_0$ and the constants $P$, $Q$, $J$ and $s$, which map onto Levi-Civita's exponents only through Gleiser's (32), (36) and (38), and its metric functions carry irrational powers of two linear factors; the page keeps $\sigma$ and the four charts above cover both signs.
The charts of single members, Griffiths and Podolský's $p$ at $\sigma = 0$, da Silva and coauthors' $X$ at $\sigma = \pm 1/2$ and Žofka and Bičák's black string coordinates at $\sigma = 1/2$, belong to the pages of those members, the anti-de Sitter soliton and the topological black holes.

## How the charts are printed

The checker's canonical form proves an identity between symbolic powers whose numeric parts are integers, and not between $(1 - \Lambda y)^{5/3 - p}$ and $(1 - \Lambda y)^{2/3 - p}$, so no symbolic exponent on the page carries a fraction.
In the two proper distance charts $\cos^{4/3}$ and $\cosh^{4/3}$ of $\sqrt{3|\Lambda|}\rho/2$ are both $(1 + 3\Lambda P^2/4)^{-2/3}$, and $P$ is held with the slope $dP/d\rho = 1 + 3\Lambda P^2/4$, $\sec^2$ or $\mathrm{sech}^2$ of half the argument; in Gleiser's chart $Y$ is held with the slope $(1 + \Lambda Y)^2$, and in Sarıoğlu and Tekin's $h = 1 - r_0^3/r^3$ with the slope $3(1 - h)/r$, which holds no $r_0$.
So every value is rational in the held name, and in $r$ for Sarıoğlu and Tekin's chart, times one power with a symbolic exponent of $P$, $Y$ or $h$.
Written out, $y^p(1 - \Lambda y)^{2/3 - p}$ fails both ways: the checker's form does not merge $(1 - \Lambda y)^{5/3 - p}$ with $(1 - \Lambda y)^{2/3 - p}$, and $(y/(1 - \Lambda y))^p$ comes back with powers of $-1$ carrying the exponents.
$4 + 3\Lambda P^2$ is printed as four times the sum $1 + 3\Lambda P^2/4$ the line element names, so that its fractional powers carry no radical constant.
`lt_powers` in `print_charts.py` prints that power whole on its name and factors the rest, holding the power as a placeholder so that the factoring never merges it back.

The tensors are printed for free exponents, which holds on the circle.
On the circle the Ricci tensor is $\Lambda g_{\mu\nu}$, the Einstein tensor $-\Lambda g_{\mu\nu}$, the Ricci scalar $4\Lambda$, and the Weyl tensor the Riemann tensor less $\Lambda(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho})/3$; `linet_tian_on_circle` checks each on the checker's parametrisation of the circle and writes them so.

## The Kretschmann scalar

On the circle, with $u = \sqrt{3\Lambda}\rho/2$,

$$K = \frac{8\Lambda^2}{3} + \frac{4\Lambda^2}{3\cos^4u} - \frac{9p_0p_2p_3\Lambda^2\cos 2u}{\sin^4u\cos^4u} = \frac{8\Lambda^2}{3} + \frac{4\Lambda^2}{3\cos^4u} - \frac{16p_0p_2p_3\cos(\sqrt{3\Lambda}\rho)}{Q^4},$$

found by fitting the computed scalar in $\sin^2u$ and $\sigma$, where the numerator came out $S^2\Sigma^3 - 54Ss^2(2s - 1)^2 + 27s^2(2s - 1)^2$ with $S = \sin^2u$ and $p_0p_2p_3 = -4s^2(2s - 1)^2/\Sigma^3$.
As $\Lambda \to 0$ it is Levi-Civita's $-16p_0p_2p_3/\rho^4$.
Continued to $\Lambda < 0$ and written in Sarıoğlu and Tekin's $h = \tanh^2u$ and $r = r_0\cosh^{2/3}u$, it is their $24/L^4 + 12r_0^6/(L^4r^6) - 81(r_0^9/(L^4r^9))((1 + h)/h^2)\prod p_k$ term by term, a check from a source.
The excess over $8\Lambda^2/3$ is the square of the Weyl tensor, for this static metric a sum of squares of the electric part: $\Lambda^2\cos^{-4}u\,(4/3 - 9p_0p_2p_3\cos 2u/\sin^4u)$ with $-4/27 \le p_0p_2p_3 \le 0$, which vanishes nowhere in the open range of $\rho$ for either sign of $\Lambda$, so no member is de Sitter or anti-de Sitter space.
The axis $\rho = 0$ is regular only where $p_0p_2p_3 = 0$, $\sigma = 0$ or $1/2$, and the far end $\rho = \pi/\sqrt{3\Lambda}$ only where $p_0p_2p_3 = -4/27$, $\sigma = -1/2$, $1/4$ or $1$, as da Silva and coauthors' (4) has it.

## Units

A power whose exponent is not a number is read with its base in a fixed unit, so $P^{2p}$, $(y/(1 - \Lambda y))^p$ and $h^p$ are pure numbers on the page; the conicities $B$ and $C$ then carry an inverse length, as Žofka and Bičák say of their $1/C$.
