# The double Kerr solution

Kramer and Neugebauer's two Kerr black holes on one axis, in the canonical chart of Weyl, Lewis, and Papapetrou, are

$$ds^2 = -f\left(c\,dt - \omega\,d\phi\right)^2 + \dfrac{1}{f}\left(e^{2\gamma}\left(d\rho^2 + dz^2\right) + \rho^2d\phi^2\right),$$

with $f$, $\omega$ and $\gamma$ functions of $\rho$ and $z$.
The chart is written by `_tools/derivations/print_charts.py --metric double_kerr`, in about a minute, and `verify_metrics.py --system double_kerr/weyl` checks it in seconds.

## Step 1. Why the three functions are left free

The literature on the solution uses this one chart; Kramer and Neugebauer, Tomimatsu and Kihara, Manko and Ruiz, Neugebauer and Hennig, and Herdeiro and Rebelo all write the metric in it.
The solution's own $f$, $\omega$ and $\gamma$ are quotients of determinants of the four distances $r_j = \sqrt{\rho^2 + (z - a_j)^2}$, and a curvature component written out in them runs to tens of thousands of terms: the static pair alone, Bach and Weyl's two Schwarzschild rods, did not finish its metric in the checker's canonical form in two minutes.
So the chart leaves the three functions free, as Gowdy's and Einstein and Rosen's charts do, no component assumes a field equation, and the parameters' descriptions state the field equations and the solution.
`OVERRULED` in `metric_tags.py` sets the tag `vacuum`, which the free functions cannot show.

## Step 2. The field equations

With $x^0 = ct$ the function $\omega$ is a length.
The vacuum equations are Ernst's equation for $\mathcal{E} = f + i\chi$, $f\,\nabla^2\mathcal{E} = (\nabla\mathcal{E})^2$ in the flat space of $\rho$, $\phi$ and $z$, with the twist potential $\partial_\rho\omega = -\rho f^{-2}\partial_z\chi$ and $\partial_z\omega = \rho f^{-2}\partial_\rho\chi$, and a quadrature for $\gamma$.
Written in $f$ and $\omega$ they are the three equations the parameters state, which do not depend on the sign of $\omega$.
`double_kerr_check` replaces $\partial_z^2f$, $\partial_z^2\omega$ and both first derivatives of $\gamma$ by them and holds every component of the Ricci tensor to vanishing, and holds the quadrature to being integrable.

## Step 3. Kramer and Neugebauer's solution

The Ernst potential is $\mathcal{E} = D_-/D_+$ with $D_\pm$ the determinant of $X_{ik} \pm 1$, $X_{ik} = (\alpha_ir_i + \alpha_kr_k)/(a_i - a_k)$, over $i = 1, 3$ and $k = 2, 4$, the form Herdeiro and Rebelo give after Yamazaki, their section 2.1, with their $e^{i\omega_j}$ written $\alpha_j$ since $\omega$ is the metric function here.
Kramer and Neugebauer's own constants are four positions $K_l$ and four phases $\omega_l$, as Manko and Ruiz's comment of 2005 records.
Then $f = \mathrm{Re}(D_-D_+^*)/|D_+|^2$, $e^{2\gamma} = \mathrm{Re}(D_-D_+^*)/(C\,r_1r_2r_3r_4)$, and $\omega = \mathrm{Im}(M^*D_+ - LD_+^*)/\mathrm{Re}(D_-D_+^*)$ with Yamazaki's potentials $L$ and $M$.
`double_kerr_check` holds the potential to Ernst's equation and the stated $e^{2\gamma}$ to the quadrature, at random rod ends and phases, in forty digits.
Inside an ergoregion $f < 0$ and $e^{2\gamma}$ is negative with it, so $g_{\rho\rho}$ stays positive and $\gamma$ has the imaginary part $\pi/2$, which the description of $\gamma$ says.

## Step 4. The pair that is drawn

Every drawing is Herdeiro and Rebelo's pair of equal masses $M$ and opposite angular momenta $\pm J$, their section 3, at $J = M^2$ and $\zeta = 4M$ in units with $G = c = 1$.
Their (3.6) gives each rod the length $a = 2\sqrt{2/3}\,M$, and the tangents of half their phases are $b = -(3 - \sqrt{6})/3$ and $c = 3 - \sqrt{6}$.
A turn of every phase by $\pi$ exchanges $D_-$ and $D_+$ and so the sign of the masses; read with $b$ and $c$ as those tangents the quotient $D_-/D_+$ has $f \to 1 + 4M/r$, so the drawings take $\alpha_j = -e^{i\omega_j}$, for which $f \to 1 - 4M/r$.
`_tools/test_double_kerr.py` holds the pair, in complex arithmetic no drawing uses, to their closed forms: Komar mass $1$ and angular momentum $\pm 1$ for each hole, $\omega = 1/\Omega_H$ on each rod with $\Omega_H = (2M - a)/4J$, the surface gravity $a(2M - a)/4M^2$ at each pole, the force $1/12$ on the strut, $\omega = 0$ on the axis and no NUT charge.
At $J = M^2$ a lone Kerr hole would be extreme; in the pair the limit is $|J| = \sqrt{3}\,M^2$, their (1.6).

## Step 5. The diagrams

The rotations fix the axis, and the half turn $z \to -z$, $\phi \to -\phi$, which exchanges the holes, fixes the half plane $z = 0$, $\phi = 0$, so both are totally geodesic and their null curves are null geodesics.
On both $\omega = 0$, and the metric is $-f\,c^2dt^2 + (e^{2\gamma}/f)$ times $dz^2$ or $d\rho^2$.
The spacetime diagrams are drawn from real expressions for the three functions, every complex number carried as two parts, and the midplane view checks them against the published Einstein tensor; `--verify` checks the rays against the tortoise coordinate taken from the complex arithmetic.
On the axis the published $g^{tt}$ is $0/0$, so that view orients its cones by $\partial_t$, and the two rods, each a whole horizon, are hatched with their four poles marked.
On the plane $z = 0$ the four distances are two, $U$ and $V$, and $f$ and $e^{2\gamma}$ are quotients of quartics in them with $C = 36/25$, which the embedding and the conformal diagram use; at $\rho = 0$, $f = 5/41$ and $e^\gamma = 3/4$.
The embedded plane is a cone of angle $\tfrac{4}{3}\cdot 2\pi$ at the strut, drawn in Minkowski space as far as $\rho = 1.09\,m$, where it lies level, and in flat space beyond.
The conformal diagram draws the axis above the holes, a diamond with the horizon's pole on its left, the strut, a diamond with a horizon on every side, and the plane $z = 0$, Minkowski's triangle.
