# Kerr-de Sitter and Kerr-anti-de Sitter

Carter's rotating black hole with a cosmological constant of either sign,

$$ds^2 = -\frac{\Delta_r}{\rho^2}\left(c\,dt - \frac{a\sin^2\theta}{\Xi}d\phi\right)^2 + \frac{\rho^2}{\Delta_r}dr^2 + \frac{\rho^2}{\Delta_\theta}d\theta^2 + \frac{\Delta_\theta\sin^2\theta}{\rho^2}\left(a\,c\,dt - \frac{r^2 + a^2}{\Xi}d\phi\right)^2,$$

$$\Delta_r = \left(r^2 + a^2\right)\left(1 - \frac{\Lambda r^2}{3}\right) - r_s r,\qquad \Delta_\theta = 1 + \frac{\Lambda a^2}{3}\cos^2\theta,\qquad \Xi = 1 + \frac{\Lambda a^2}{3},\qquad \rho^2 = r^2 + a^2\cos^2\theta.$$

Its five charts are written by `_tools/derivations/print_charts.py --metric kerr_de_sitter`, which took 229 seconds on 1 October 2026, and `verify_metrics.py --system kerr_de_sitter/<chart>` checks each in under a minute.

## Step 1. Where the metric comes from

Carter's paper of 1968 (`carter1968separable`) gives the solution as his family [A], equations (4) to (6), in polynomial coordinates and with the opposite sign of $\Lambda$; it has no $\Xi$ and never writes this member out.
The form above is equation (4.1) of Hawking, Hunter and Taylor-Robinson (`hawking1999rotation`) with their $l^2 = -\Lambda/3$, and equation (1) of Caldarelli, Cognola and Klemm and (2.1) of Gibbons, Perry and Pope: only $d\phi$ is divided by $\Xi$.
Akcay and Matzner, Stuchlík and Slaný, Lake and Zannias, Borthwick, and Dias, Eperon, Reall and Santos divide $dt$ by $\Xi$ as well, so their $t$ is $\Xi$ times this one, and every formula taken from them here is rescaled.
The mass enters as the length $r_s = 2GM/c^2$, as Kottler's does on the site, so that $a = 0$ is `schwarzschild_de_sitter` term by term and $\Lambda = 0$ is `kerr` with $2GM/c^2$ written $r_s$.

## Step 2. The four names

Every chart declares $\Xi$, $\rho$, $\Delta_r$ and $\Delta_\theta$ as names for the expressions above, the way Weyl's chart of the Curzon-Chazy particle declares $R$.
`Reader` reads a name spelled with a command in its subscript, `\Delta_\theta`, since 1 October 2026, and knows `\Xi`.
The checker hands every value back in $r$, $a$, $\Lambda$, $r_s$, $\sin\theta$ and $\cos\theta$, and `KerrDeSitterForms` in `print_charts.py` writes it in the names: it factors the value with the even powers of the sine written in the cosine, writes each factor that is a multiple of a named polynomial by its name, writes $(u + v)(u - v)$ as $u^2 - v^2$, and writes every sum left over in the shortest of five forms, as it stands in the cosine or the sine, grouped by the powers of $r_s$, of $\Lambda$ or of both, or grouped by the powers of $\Delta_r$ after $r_s r$ is written as $(r^2 + a^2)(1 - \Lambda r^2/3) - \Delta_r$.
That is how $g_{tt}$ comes out as $-(\Delta_r - a^2\Delta_\theta\sin^2\theta)/\rho^2$ and $\Gamma^t{}_{tr}$ with the numerator $3r_s(r^2 - a^2\cos^2\theta) - 2\Lambda r\rho^4$.
Every printed value is read back and compared with the value it came from, so the choice of form cannot change a value.

## Step 3. The curvature

The spacetime is an Einstein space, $R_{\mu\nu} = \Lambda g_{\mu\nu}$ and $R = 4\Lambda$, which `kerr_de_sitter_check` holds every chart to before it is written.
So the Riemann tensor is the Weyl tensor plus the constant curvature part,

$$R_{\mu\nu\rho\sigma} = C_{\mu\nu\rho\sigma} + \frac{\Lambda}{3}\left(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}\right),$$

which the same check holds in all 256 slots, and `kerr_de_sitter_riemann` writes each Riemann component as those two parts over one denominator: written out whole, $R_{trtr}$ runs to five lines.
The Weyl tensor carries no $\Lambda$ beyond the names: its one independent scalar is Kerr's, $\Psi_2 = -GM/(c^2(r - ia\cos\theta)^3)$, so
$$K = \frac{12r_s^2\left(r^2 - a^2\cos^2\theta\right)\left(\rho^4 - 16a^2r^2\cos^2\theta\right)}{\rho^{12}} + \frac{8\Lambda^2}{3},$$
Kerr's Kretschmann scalar plus de Sitter's, which the script checks against sympy.

## Step 4. The other four charts

Each is checked, slot by slot, to be the first chart pulled back through the Jacobian `kerr_de_sitter_jacobian` writes.

The chart that does not rotate at infinity keeps $t$, $r$ and $\theta$ and takes $\Phi = \phi - a\Lambda ct/3$, the angle of Henneaux and Teitelboim's (B.2a) (`henneaux1985`) and of Gibbons, Perry and Pope's (2.6) to (2.8): its $g_{t\Phi} = -ar_sr\Delta_\theta\sin^2\theta/(\Xi^2\rho^2)$ vanishes with the mass and falls off at infinity.
With $\Xi - \Lambda a^2\sin^2\theta/3 = \Delta_\theta$ and $\Xi - \Lambda(r^2 + a^2)/3 = 1 - \Lambda r^2/3$ the two one forms become $(\Delta_\theta\,c\,dt - a\sin^2\theta\,d\Phi)/\Xi$ and $(a(1 - \Lambda r^2/3)\,c\,dt - (r^2 + a^2)\,d\Phi)/\Xi$.

The Kerr charts follow the principal null congruences, $dv = c\,dt + (r^2 + a^2)\,dr/\Delta_r$ and $d\tilde\phi = d\phi + a\Xi\,dr/\Delta_r$ for the ingoing one and the opposite signs for the outgoing one.
They are Borthwick's (50) and (51) (`borthwick2018`) and Dias, Eperon, Reall and Santos's (3.9) and (3.10) (`dias2018`), whose time is $\Xi$ times this one.
Then $c\,dt - a\sin^2\theta\,d\phi/\Xi = (dv - a\sin^2\theta\,d\tilde\phi/\Xi) - \rho^2dr/\Delta_r$ and $a\,c\,dt - (r^2 + a^2)\,d\phi/\Xi = a\,dv - (r^2 + a^2)\,d\tilde\phi/\Xi$, so $g_{rr} = 0$ and the chart runs through every root of $\Delta_r$.
The diagrams fix the constants by $v = ct + r_*$ and $u = ct - r_*$ with $r_* = 0$ at $r = 0$.

The Kerr-Schild chart is Appendix A, (A.1) to (A.3), of Gibbons, Lü, Page and Pope (`gibbons2005kds`) with their $\lambda = \Lambda/3$ and $2M = r_s$: the de Sitter metric in spheroidal coordinates plus $r_sr/\rho^2$ times the square of a null one form, reached by their (3.3), $c\,d\tau = c\,dt + r_sr\,dr/((1 - \Lambda r^2/3)\Delta_r)$ and $d\psi = d\phi - (a\Lambda/3)\,c\,dt + ar_sr\,dr/((r^2 + a^2)\Delta_r)$.
Its time and angle are those of the frame that does not rotate at infinity, so for $\Lambda > 0$ it ends at $r = \sqrt{3/\Lambda}$, where the static chart of the de Sitter background does, beyond $r_c$.
Its metric has three cross terms, and sympy's `adjugate()` and `det()` did not finish on it in ten minutes; `_adjugate` in `verify_metrics.py` takes both by cofactors with every minor put through `norm`, in under a second, and the whole collection agrees with sympy as before.

## Step 5. The values the diagrams are drawn at

With $\Lambda > 0$: $r_s = 1$, $a = 0.45\,r_s$, Kerr's $0.9\,GM/c^2$, and $\Lambda = 0.2/r_s^2$, Kottler's.
$\Delta_r$ has the roots $r_n = -4.2958$, $r_- = 0.2788$, $r_+ = 0.7848$ and $r_c = 3.2323$, with surface gravities $\kappa = |\Delta_r'|/(2(r^2 + a^2))$ of $0.313$, $0.813$, $0.256$ and $0.170$ on the axis, and on the equator $g_{tt} = 0$ where $\Delta_r = a^2$, at $1.1048$ and $3.1734$.
With $\Lambda < 0$: $\ell = \sqrt{-3/\Lambda} = 1$, $r_s = 2\ell$, Hawking and Page's, and $a = \ell/2$, so that $r_- = 0.1369$, $r_+ = 0.8594$ and the ergosurface on the equator is at $0.9386$.
At $a = 0.6\,\ell$ the circle of least circumference on the equator, $r^3 = r_sa^2/2$, lies outside $r_+$, and the embedded surface would narrow beyond its throat; at $a = \ell/2$ it lies inside.

On the axis the metric is $-f\,c^2dt^2 + dr^2/f$ with $f = \Delta_r/(r^2 + a^2)$, and $1/f$ has no polynomial part, so $r_* = \sum_i \ln|1 - r/r_i|/f'(r_i)$ over the four roots, real or complex, which vanishes at $r = 0$.
For $\Lambda > 0$ it tends to the same $R$ at $r \to +\infty$ and $r \to -\infty$, since the residues sum to zero; for $\Lambda < 0$ the complex pair gives two different limits.
`CarterAxis` and `CarterAdSTower` in `conformal.py` carry the maps, and their docstrings the derivation.
