# The naked black holes of Horowitz and Ross: derivation notes

The page holds three charts, each from Gary T. Horowitz and Simon F. Ross, "Naked black holes", Phys. Rev. D 56, 2180 (1997), hep-th/9704058 (HR97), and "Properties of naked black holes", Phys. Rev. D 57, 1098 (1998), hep-th/9709050 (HR98), read from the arXiv TeX sources.
Both articles work in Planck units for the general relativity examples and in string units, $\alpha' = 1$, for the string example; the page reads every parameter and $r$ as a length and keeps $c$ in the line element, with $g_{tt}$ dimensionless.
`naked_black_holes.py` beside this file holds the charts and their checks, and `naked_black_holes_drawings.py` the drawings.

## The charts and their sources

- `dilaton`, "Dilaton Coupling $a$": HR97 eqs. (gensol), (uuf), and (R), $ds^2 = -F\,dt^2 + dr^2/F + R^2d\Omega^2$ with $F = (r - r_+)(r - r_-)/R^2$ and $R = r(1 - r_-/r)^{a^2/(1+a^2)}$, the dilaton black holes of Gibbons and Maeda and of Garfinkle, Horowitz, and Strominger for the action (dbhac), $R - 2(\nabla\phi)^2 - e^{-2a\phi}F^2$.
  The horizon is $r_+$, and $r_-$ is a curvature singularity for $a \neq 0$; HR97 say so beside (R).
  The magnetic hole has $Q^2 = r_+r_-/(1 + a^2)$, HR97 (dbhmass), and $e^{-2\phi} = (1 - r_-/r)^{2a/(1+a^2)}$, HR98 (dbhdil), whose TeX source carries a comment asking for a check; the field equations hold at points, so the page states it.
  The chart holds $R$ as a function of $r$ with its slope $R' = R((1 + a^2)r - r_-)/((1 + a^2)r(r - r_-))$, since sympy does not finish factoring powers in $1/(1 + a^2)$, as the Chan-Horne-Mann chart of `dilaton_unusual_asymptotics` found first; the Ricci scalar prints with no "R =" before it, since $R$ is the areal radius here.
  The domain takes $a > 0$, so that $r_-$ is a singularity, and the parameter description says that $a = 0$ is Reissner and Nordström's hole with an inner horizon at $r_-$.
- `two_charges`, "Electric and Magnetic Charges": HR97 eqs. (uf), (ur), and (mcff), the same $F$ with $R^2 = r^2 - \Sigma^2$, for the action (uact), $R - 2(\nabla\phi)^2 - e^{-2\phi}(F^2 + G^2)$, with electric charge $Q$ in $F$ and magnetic charge $P$ in $G$, and with $\Sigma = (P^2 - Q^2)/2M$ and $r_\pm = M \pm (M^2 + \Sigma^2 - Q^2 - P^2)^{1/2}$.
  The field equations hold with $e^{2\phi} = (r + \Sigma)/(r - \Sigma)$, $F_{tr} = -Qe^{2\phi}/R^2$, and $G = P\sin\theta\,d\theta\wedge d\phi$, and fail with the opposite sign of $\phi$.
  The domain is $r > |\Sigma|$ and $|\Sigma| \le r_- < r_+$: $r_+$ is the event horizon, $r_-$ the inner horizon, and $|\Sigma|$ the singularity, all three named by HR97 beside (mcff).
- `string_ns`, "String Metric, Winding and Momentum": HR97 eqs. (nsmet), (dd), and (nsdil), after Sen, $ds^2 = -\Delta^{-1}(1 - r_0/r)dt^2 + (1 - r_0/r)^{-1}dr^2 + r^2d\Omega^2$ with $\Delta = (1 + r_0\sinh^2\gamma_1/r)(1 + r_0\sinh^2\gamma_p/r)$ and $e^{2\phi} = \Delta^{-1/2}$.
  The page writes it in the lengths $r_1 = r_0\sinh^2\gamma_1$ and $r_p = r_0\sinh^2\gamma_p$, so $\Delta = (r + r_1)(r + r_p)/r^2$, $-g_{tt} = r(r - r_0)/((r + r_1)(r + r_p))$, and every value is rational in $r$; printed in $\gamma_1$ and $\gamma_p$, sympy wrote each $\sinh$ in exponentials and the chart ran to 390,000 characters.
  The charges grow as $\sinh 2\gamma_1$ and $\sinh 2\gamma_p$, HR97 (nscharge), which the parameter descriptions say.
  The Ricci scalar near $r = 0$ grows as $r_0/r^3$ for $r_1, r_p > 0$, so $r = 0$ is a curvature singularity, as the domain says.

No chart is a repeat of a live page.
At $a = 1$ the dilaton chart is the static chart of `dilaton_black_hole` with $r_s = r_+$ and $r_d = r_-$; at $a = \sqrt{3}$ it is the chart "Einstein Metric, Four Dimensions" of `kaluza_klein_black_hole` with its $r$ the page's $r - r_-$, $r_s = r_+ - r_-$, and $q = r_+$; at $a = 0$ it is `rn_metric` with $r_s = r_+ + r_-$ and $r_q^2 = r_+r_-$; and at $r_- = 0$ it is `schwarzschild`.
The chart with two charges is `rn_metric` at $\Sigma = 0$, $Q = P$, and the $a = 1$ hole with its $r$ the page's $r + r_-$, $r_s = r_+ + r_-$, and $r_d = 2r_-$, at $|\Sigma| = r_-$, one charge zero, as HR97 say beside (mcff).
The Einstein metric $e^{-2\phi}g = \Delta^{1/2}g$ of the string chart is the $a = 1$ hole for $r_1 = r_p$, with its $r$ the page's $r + r_1$, $r_s = r_0 + r_1 = r_0\cosh^2\gamma$, and $r_d = r_1$, and the $a = \sqrt{3}$ hole of `kaluza_klein_black_hole` for $r_p = 0$, with $r_s = r_0$ and $q = r_0 + r_1$, as HR97 say at the end of their Sec. 4.1.
Each of these is checked at points against the live chart's own line element.

## The tidal components

For $ds^2 = -F\,dt^2 + dr^2/F + R^2d\Omega^2$ the static orthonormal frame has $R_{\hat0\hat2\hat0\hat2} = F'R'/(2R)$ everywhere and $R_{\hat2\hat3\hat2\hat3} = (1 - FR'^2)/R^2$, which is HR97's (stcur), $1/R^2$, at the horizon, where $F = 0$.
A radial fall of energy $E$ per unit mass has $u^t = E/F$ and $u^r = -(E^2 - F)^{1/2}$, a radial boost of the static frame with $\cosh\alpha = E/F^{1/2}$, HR97 (inframea), and the transverse tidal component in the falling frame is HR97's (bhtide), $R_{0'20'2} = -[R''(E^2 - F) - F'R'/2]/R$.
HR97's (gdevb), with $g_{tt} = -F/G$, gives for the string chart $R_{0'20'2} = -(E^2G' - F')/(2r)$ with $F = 1 - r_0/r$, $G = \Delta$, and $R = r$, and at $r_0$ their (nsgeod), $-E^2\Delta'(r_0)/(2r_0) + 1/(2r_0^2)$, with $\Delta'(r_0)$ their (dd).
The boost law, $R_{0'k0'k} = R_{0k0k} + \sinh^2\alpha(R_{0k0k} + R_{1k1k})$, HR97 (bcurv), and both components were derived again in sympy for the task with $F$, $G$, and $R$ free functions of $r$ before the History went up.
The chart checks compute the frame components from each chart's own Riemann tensor at points and hold them to (stcur) and (bhtide), and the string chart's to (gdevb), (nsgeod), and (dd).
At the horizon $R''/R$ is $-a^2(1 - \epsilon)^2/((1 + a^2)^2r_+^2\epsilon^2)$ with $\epsilon = 1 - r_-/r_+$, HR97 (dbhtide), and $-(1 - \delta)/(r_+^2\delta^2)$ with $\delta = 1 - \Sigma^2/r_+^2$, HR97 (utide); HR97 print the magnitudes, and the checks hold both.

## The checks

`chart_check` in `naked_black_holes.py` runs in `print_charts.py` before anything is written, in forty digits at points chosen off every special value, with no symbolic simplification.

- Each chart's metric against its source equations.
- The dilaton chart: Einstein's, Maxwell's, and the dilaton's equations for the magnetic hole of the action (dbhac), and the mass $M = r_+/2 + (1 - a^2)r_-/(2(1 + a^2))$ of HR97 (dbhmass) from $r(1 + g_{tt})$ at $r = 10^{30}$.
- The chart with two charges: the same equations for the action (uact) with both fields, the mass $M = (r_+ + r_-)/2$, and $r_\pm$ of (mcff) from the charges.
- The tidal components and the special cases above.

The three charts made 85, 92, and 22 comparisons on 2026-11-23, the largest residual $2 \times 10^{-41}$.

## The drawings

Every drawing takes one member of each chart, in units of its horizon, and `naked_black_holes_drawings.py` records each choice in its docstring.

- The dilaton chart is drawn at $a = 1/2$ and $r_- = r_+/2$, a coupling no live page draws: $a = 1$ is `dilaton_black_hole` and $a = \sqrt{3}$ is `kaluza_klein_black_hole`.
  Its tortoise coordinate $r_* = \int dr/F$ has no closed form there, so it is summed by quadrature, split as $A\ln|(r - r_+)/(r_+ - r_-)|$ with $A = R(r_+)^2/(r_+ - r_-)$ plus the integral of a function smooth at $r_+$, taken from the singularity, where the substitution $r - r_- = w^{5/2}$ removes the integrable power $(r - r_-)^{-3/5}$.
  The sum agrees with $1/F$ to $4 \times 10^{-9}$ by differences, and $r_* = 0$ at the singularity puts it on $UV = 1$ of Kruskal's coordinates with the surface gravity $\kappa = 1/(2A)$.
- The chart with two charges is drawn at $r_- = 0.6\,r_+$ and $\Sigma = 0.3\,r_+$, so $\kappa_-/\kappa_+ = (r_+^2 - \Sigma^2)/(r_-^2 - \Sigma^2) = 91/27$, and its conformal diagram is the tower of `rn_metric`, $r_*$ shifted by a constant to vanish at the singularity $r = |\Sigma|$, which puts the singularity on the vertical lines $X = \pm\pi/2$ that the tower draws; the label of the singularity is written $r = |\Sigma|$.
- The string chart is drawn at $r_1 = r_p = r_0/2$, where $r_* = r + (r_0 + r_1)\ln|r/r_0 - 1|$ in closed form, zero at the singularity $r = 0$, with $\kappa = 1/(2(r_0 + r_1))$.
- Each chart has a spacetime diagram of $t$ and $r$, oriented by the ingoing rays, as the static charts of `dilaton_black_hole` and `rn_metric` are.
- Each chart has a conformal diagram: Schwarzschild's square for the dilaton and string charts, with a spacelike singularity, and the tower for the chart with two charges, with a timelike one.
- The embedding diagram draws the equatorial plane of the dilaton chart's slice $t = 0$ at the member drawn, its circles of radius $R$, from the throat at the bifurcation sphere $r = r_+$ out to $6r_+$ on both exteriors.
  The other two charts are left out of it: the chart with two charges is another black hole than the one embedded, and the string chart measures the manifold with the string metric; the slices table says so for each.
  The moment $t = 0$ is drawn on the dilaton chart's spacetime diagram and conformal diagram.
- The curvature figure draws the stored Kretschmann and Ricci scalars of the dilaton chart at the member drawn against $r$, from $0.55\,r_+$ to $6r_+$, a box of two panels.
  The figure writes the Ricci scalar $\mathcal{R}$, since $R$ is the areal radius on this page.
  Its checks hold the Ricci scalar to $2(\nabla\phi)^2$, the trace of Einstein's equations with the dilaton of HR98 (dbhdil), the Maxwell stress being traceless.
- No view is a movie, so the page adds nothing to the MOVIES or STACKS tables.

Nothing is left undrawn.
