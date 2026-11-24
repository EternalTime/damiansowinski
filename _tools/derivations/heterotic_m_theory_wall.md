# The domain wall of heterotic M-theory: derivation notes

The page holds four charts of the walls at the ends of the interval of Hořava and Witten's theory, reduced on a Calabi-Yau space to gauged $N = 1$ supergravity in five dimensions by André Lukas, Burt A. Ovrut, Kellogg S. Stelle, and Daniel Waldram, "Universe as a domain wall", Phys. Rev. D 59, 086001 (1999), hep-th/9803235 (LOSW), read from the arXiv TeX source.
The scalar of the solution is the volume $V$ of the Calabi-Yau space, a pure number in units of a fixed reference volume, LOSW's eq. (metric1), and $\phi = \ln V$.
Every chart is in the signature $(-,+,+,+,+)$, $\alpha$ is LOSW's constant, an inverse length, and $c$ stands in the line element alone.
`heterotic_m_theory_wall.py` beside this file holds the charts and their checks, and `heterotic_m_theory_wall_drawings.py` the drawings.

## The charts and their sources

- `static`, "Lukas, Ovrut, Stelle, and Waldram, Static Walls": LOSW eq. (sol), $ds^2 = a_0^2H\eta_{\mu\nu}dx^\mu dx^\nu + b_0^2H^4dy^2$, $H = \frac{\sqrt{2}}{3}\alpha|y| + c_0$, $V = b_0H^3$, on $-\pi\rho \le y \le \pi\rho$ with $\pi\rho$ and $-\pi\rho$ one point and $y \to -y$ a symmetry.
  LOSW write $\partial_y^2H = \frac{2\sqrt{2}}{3}\alpha(\delta(y) - \delta(y - \pi\rho))$, so the wall at $y = 0$ is the delta of $|y|$ and the wall at $y = \pm\pi\rho$ has the opposite delta; the components are those of $|y|$, which hold for $|y| < \pi\rho$.
- `conformal`, "Lehners, McFadden, and Turok, Conformally Flat": Jean-Luc Lehners, Paul McFadden, and Neil Turok, Phys. Rev. D 75, 103510 (2007), hep-th/0611259 (LMT), eq. (domainwall), $ds^2 = h^{2/5}(\eta_{\mu\nu}dx^\mu dx^\nu + dy^2)$, $e^\phi = h^{6/5}$, $h = 5\alpha_{\rm LMT}(y + y_0) + c$ on $-y_0 \le y \le y_0$.
  Their footnote gives $\alpha_{\rm LMT} = \alpha/3\sqrt{2}$, so the slope is written $\frac{5\sqrt{2}}{6}\alpha$; their $c$, which on this page is the speed of light, is written $h_0$, the value of $h$ on the wall at $w = -w_0$, and their $y$ and $y_0$ are written $w$ and $w_0$, so that one letter means one thing.
  LMT put the wall of negative tension at $-y_0$ and the wall of positive tension at $+y_0$.
  With $b_0 = a_0^6$ the static chart is this one through $dw/dy = (b_0/a_0)H^{3/2}$ and $h = a_0^5H^{5/2}$, which the checks hold.
- `static_bulk`, "Chamblin and Reall, Static Bulk": the bulk found by H. A. Chamblin and H. S. Reall, Nucl. Phys. B 562, 133 (1999), hep-th/9903225 (CR), Sec. 8, in LMT's form, their eq. (staticbulk), $ds^2 = -f\,dt_{\rm LMT}^2 + r^{12}dr^2/f + r^2d\vec{x}^2$, $f = \alpha_{\rm LMT}^2r^2 - \mu r^4$, $e^\phi = r^6$.
  Their time has the dimensions of a length squared, so it is written $t_{\rm LMT} = (c/\alpha_{\rm LMT})t$ with $t$ a time, and $\mu$ is written $\lambda^2\alpha_{\rm LMT}^2$ with their own $\lambda = \sqrt{\mu}/\alpha_{\rm LMT}$ of Appendix C: $ds^2 = -r^2(1 - \lambda^2r^2)c^2dt^2 + 18r^{10}dr^2/\alpha^2(1 - \lambda^2r^2) + r^2d\vec{x}^2$.
  LMT print the range as $0 \le r \le \sqrt{\mu}/\alpha$, but their $f$ vanishes at $r = \alpha/\sqrt{\mu}$, and their own Appendix C writes the static region as $\bar r = \lambda r \le 1$; so the domain is $0 < r < 1/\lambda$, with the timelike naked singularity at $r = 0$ and the cosmological horizon, CR's name for it, at $r = 1/\lambda$.
  LMT recover the static wall at $\mu = 0$, and CR's walls stay at rest wherever they are placed at $M = 0$; CR's $M$ and LMT's $\mu$ differ by a positive factor, so $\lambda = 0$ is $M = 0$.
- `rolling`, "Lukas, Ovrut, and Waldram, Rolling Walls": Lukas, Ovrut, and Waldram, Phys. Rev. D 60, 086001 (1999), hep-th/9806022 (LOW), eqs. (burt6), (burt15), and (burt16) in the gauge $n = 1$ of (burt14): $N = a(y)$, the walls scaled by $A(t - t_0)^p$ and the interval and $V$ by $B(t - t_0)^q$, $p = \frac{3}{11}(1 \mp \frac{4}{3\sqrt{3}})$, $q = \frac{2}{11}(1 \pm 2\sqrt{3})$, with LOSW's profile in $y$.
  The chart covers the branch $t > t_0$, and the scale factors are held as names, $a = A(t - t_0)^p$ and $b = B(t - t_0)^q$, with declared rates $\dot a = pa/(t - t_0)$ and $\dot b = qb/(t - t_0)$, so every value is rational in $a$, $b$, $t - t_0$, $H$, and the parameters, as Newman and Tamburino's chart holds its names; the components are built in $x^0 = ct$, where the names are written in $x^0/c$.
  The printer's fallback for a bracketed sum is handed the value in $x^0$ itself, so this chart writes the time as $ct$ before it brackets a sum.

Every value of the two charts in $y$ is printed in $H$, $\mathrm{sgn}(y)$, and $\delta(y)$, with the coefficient of the delta taken on the wall, where $H = c_0$, as Randall and Sundrum's charts take theirs, and the Kretschmann scalar is stated for $y \neq 0$.

## The checks

Each chart is checked at points in thirty digits, apart from the printed tensors, in `heterotic_m_theory_wall.chart_check`.
The bulk is LMT's action (Action5d), $R - \frac{1}{2}(\partial\phi)^2 - 6\alpha_{\rm LMT}^2e^{-2\phi}$, which is LOSW's (actparts) with the other fields zero and their sign of the curvature turned to this one, so $R_{ab} = \frac{1}{2}\partial_a\phi\,\partial_b\phi + \frac{\alpha^2}{9}e^{-2\phi}g_{ab}$ and $\Box\phi = -\frac{2\alpha^2}{3}e^{-2\phi}$.

- Every chart solves both on both sides of the wall at $y = 0$ where it crosses it, at two sets of constants, one with $\alpha < 0$.
- On the wall at $y = 0$ the delta in $G_{\mu\nu}$ is $\sqrt{2}\alpha V^{-1}g_{\mu\nu}/\sqrt{g_{yy}}$ and in $\Box\phi$ it is $2\sqrt{2}\alpha V^{-1}/\sqrt{g_{yy}}$, LMT's boundary terms (eom1) and (eom2), $6\alpha_{\rm LMT}$ and $12\alpha_{\rm LMT}$; for the moving walls this holds at every time because the interval and $V$ scale alike, LOW's condition (burt8).
- The moving walls solve the bulk at both of LOW's branches and fail it at $p = 1/3$, so the check can fail.
- The conformal chart is the static chart with $b_0 = a_0^6$ pulled back along $w = 2(a_0^5H^{5/2} - h_0)/5k - w_0$, $k = \sqrt{2}\alpha/3$, and the static bulk at $\lambda = 0$ is the static chart with $b_0 = a_0^6$ pulled back along $r = a_0H^{1/2}$.

`verify_metrics.py --metric heterotic_m_theory_wall` checks every printed value of the four charts against sympy and the dimensions of every term; it took 96 seconds on 24 November 2026, and the rolling chart printed in 150 seconds.
`_tools/test_heterotic_m_theory_wall.py` holds each line element at points to its source equation, the horizon to the root of LMT's $f$, the History's exponents to LOW's numbers, and the drawings to one member and to the quoted branch.

## The drawings

The static walls are drawn at $\alpha = a_0 = b_0 = c_0 = 1$ with $\pi\rho = 1.2/\alpha$, so $H$ runs from 1 to $1 + 2\sqrt{2}/5 = 1.566$, and the conformal chart at the same member, $h_0 = 1$ and $2w_0 = W(\pi\rho)$, $W(y) = \mathrm{sgn}(y)(2/5k)(H^{5/2} - 1)$, so $w_0 = (3/5\sqrt{2})((1 + 2\sqrt{2}/5)^{5/2} - 1) = 0.877/\alpha$.
The walls are named by their place, at $y = 0$ and at $y = \pm\pi\rho$, and in the conformal chart by LMT's tensions.
LOSW put the gauge field $E_6$ on the wall at $y = 0$ and $E_8$ on the wall at $\pi\rho$, while the History has $V$ shrinking toward the hidden wall, which with $H$ linear in $|y|$ places the hidden wall at $y = \pi\rho$ only for $\alpha < 0$; the drawings are at $\alpha = 1$, so no wall of a drawing is named by its gauge field.
The window $\pi\rho = 1.2/\alpha$ keeps the embedding diagram of both copies of the interval inside 2:1.

- Spacetime diagrams: the plane of $t$ and the coordinate across the interval for each chart, $x_1 = x_2 = x_3 = 0$, no Christoffel symbol turning a ray out of it.
  The rays keep $ct \mp W(y)$ in the static chart, $ct \mp w$ in the conformal chart, $ct \mp r_*$ in the static bulk, $r_* = 3\sqrt{2}(32\,\mathrm{artanh}(r/2) - 16r - 4r^3/3)$ at $\lambda = 1/2$, and $\eta \mp W(y)$ for the moving walls, $\eta = (ct)^{1-q}/(1 - q)$; `null_rays.py --verify` measured a drift of at most $7 \times 10^{-7}$ along the traced rays.
  The static and moving charts draw $y$ periodic, the edges $\pm\pi\rho$ one line.
  The moving walls are drawn on the branch $p = \frac{3}{11}(1 + \frac{4}{3\sqrt{3}}) = 0.483$, $q = \frac{2}{11}(1 - 2\sqrt{3}) = -0.448$, on which the walls grow after $t_0$ while the interval shrinks, the branch LOW call closest to the early universe, at $A = B = 1$ and $t_0 = 0$.
  The static bulk's diagram ends on the horizon; the paths of the walls through it are not drawn, since LMT give them as $t(r)$, their eq. (rttraj), and no closed form of $r(t)$.
- Conformal diagrams: Minkowski's maps $p, q = \arctan(T \mp X)$ of a pair in which each plane is conformally flat.
  The static walls fill the strip $|W| \le W(\pi\rho)$ from $i^-$ to $i^+$, the conformal chart its copy $w = W - w_0$ of one side, and the moving walls, $a_0^2Hb^2(-d\eta^2 + dW^2)$ with $\eta$ finite at the singularity, the upper half of that strip with the singularity on $T = 0$.
  The static bulk, $r^2(1 - \lambda^2r^2)(-c^2dt^2 + dr_*^2)$, fills the right half of the diamond, the singularity on its axis; along a ray the affine parameter is $\int\sqrt{-g_{tt}g_{rr}}\,dr = 3\sqrt{2}\int r^6dr$, finite, $384\sqrt{2}/7$ from the singularity to the horizon at $\lambda = 1/2$, so the right edges are the horizon and no infinity.
  The region beyond the horizon, where $r$ is a time, which LMT's Appendix C extends by Kruskal's procedure, is not drawn.
- Embedding diagram: the surface of $y$ and $x_1$ at $t = 0$, $x_2 = x_3 = 0$, of the static walls, $x_1$ rolled up with period $2\pi/\alpha$ so each circle has radius $a_0H^{1/2}$, as Randall and Sundrum's strip is rolled up; $|d\rho/dl| = k/2H^{5/2} \le \sqrt{2}/6$, so the whole of it lies in flat space, two copies of the interval meeting on the wall at $y = 0$, where the slope of the profile changes sign.
  The static bulk and the moving walls are other spacetimes than the one embedded, and their drawings show no moment of it, which `exceptions/heterotic_m_theory_wall.json` records.
- Curvature figure: $K$ and the Ricci scalar against $w$ across the interval from the conformal chart at the member drawn, $K = 106\alpha^4/81h^{24/5}$ and $\mathcal{R} = 14\alpha^2/9h^{12/5}$, a plain box of two logarithmic panels.
  The static chart's Ricci scalar holds the delta of the wall at $y = 0$, which no sampled curve can draw, and the conformal chart's two scalars hold no delta; the figure checks the static chart's $K$ against it on both walls.

Every view's text is within the soft 400 characters of the shared count in `view_text.py`, the longest 383.
