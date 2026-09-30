# Majumdar-Papapetrou, extremal charged black holes at rest

Majumdar's and Papapetrou's metric is

$$ds^2 = -\frac{c^2dt^2}{U^2} + U^2\left(dx^2 + dy^2 + dz^2\right),\qquad U = 1 + \sum_i \frac{m_i}{|\mathbf{x} - \mathbf{x}_i|},$$

with $m_i = GM_i/c^2$, and every charge fixed by $GQ_i^2/(4\pi\epsilon_0c^4) = m_i^2$ and of one sign.
Its three charts are written by `_tools/derivations/print_charts.py --metric majumdar_papapetrou` in about 13 seconds, and `verify_metrics.py --system majumdar_papapetrou/<chart>` checks each.

## Step 1. The charts

The Cartesian chart is the one Majumdar, Papapetrou and Hartle and Hawking write, with $U = U(x,y,z)$ left free.
The cylindrical chart, $U = U(\rho,z)$, holds any number of holes strung along one axis, the setting of every two-centre calculation.
The isotropic chart is one hole alone, $U = 1 + m/r$, which is the extremal Reissner-Nordström black hole with areal radius $R = r + m$; its tensors are explicit, and its Kretschmann scalar is

$$K = \frac{8m^2\left(6r^2 + m^2\right)}{(r + m)^8},$$

finite at the horizon $r = 0$, where it is $8/m^4$.
The printer checks, before anything is written, that the Cartesian metric pulled back to cylindrical coordinates is the cylindrical one for any $U$, and with $U = 1 + m/\sqrt{x^2 + y^2 + z^2}$ pulled back to spherical coordinates is the isotropic one.

## Step 2. The field equations

The tensors are printed for any $U$, so no component assumes a field equation.
With $A = U^{-1}\,dt$ in units where $G = c = 4\pi\epsilon_0 = 1$, the printer checks that the Einstein tensor equals $2\left(F_{\mu\alpha}F_\nu{}^\alpha - \tfrac14 g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right)$ once $\partial_z^2U$ is replaced by $-\partial_x^2U - \partial_y^2U$, and that Maxwell's equations reduce to

$$\nabla_\mu F^{\mu t} = \frac{\partial_x^2U + \partial_y^2U + \partial_z^2U}{U^2} = 0,$$

so the Einstein-Maxwell equations hold exactly where $U$ is harmonic.
The Ricci scalar $R = -2\nabla^2U/U^3$ vanishes there too, as it must for a Maxwell field.

## Step 3. Printing

The general charts are printed with `lead = [U]` and each numerator collected by powers of $U$, and with `flip = False`, so that every value keeps the order $U\,\partial\partial U - k\,\partial U\,\partial U$ in which the Riemann components come, and $G_{xx}$ reads $-\left((\partial_xU)^2 - (\partial_yU)^2 - (\partial_zU)^2\right)/U^2$ beside $G_{yy} = \left((\partial_xU)^2 - (\partial_yU)^2 + (\partial_zU)^2\right)/U^2$.
The isotropic chart writes its metric and inverse as the line element writes $1 + m/r$, and every other value in $r + m$, with $m$ first in a product.

## Step 4. The drawings

Every two-hole drawing takes two holes of mass parameter $m$ at $z = \pm 2m$, $U = 1 + m/|\mathbf{x} - 2m\hat z| + m/|\mathbf{x} + 2m\hat z|$.
Light on the axis stays on it and light in the midplane $z = 0$ stays there, since $\partial_xU$ and $\partial_yU$ vanish on the axis and $\partial_zU$ on the midplane, so the null curves drawn on those planes are null geodesics.
On the axis $dz/dt = \pm c/U^2$, whose integral `_mp_axis` in `null_rays.py` writes in closed form in each of the three intervals the holes cut the axis into, and in the midplane

$$\int U^2\,dx = x + 4m\,\mathrm{arcsinh}\frac{x}{2m} + 2m\arctan\frac{x}{2m}.$$

The embedding diagram draws two moments $t = 0$.
One hole's equator has $dz/dr = \sqrt{2r/m + 1}\,m/r$, so

$$z = 2m\,w + m\ln\frac{w - 1}{w + 1},\qquad w = \sqrt{2r/m + 1},$$

which falls as $m\ln r$ toward the horizon while the circumference radius $r + m$ closes on $m$: the throat is infinitely long, and is drawn from $r = m/50$.
The two holes' midplane is a surface of revolution with circumference radius $\rho U$ and $dz/d\rho = \sqrt{U^2 - (U + \rho U')^2}$, real since $U + \rho U' = 1 + 8m^3/(\rho^2 + 4m^2)^{3/2}$ lies between $0$ and $U$.

The conformal diagram is one hole's maximal extension, extremal Reissner-Nordström's, drawn by `majumdar_papapetrou` in `conformal.py`.
$f = (1 - m/R)^2$ has a double root, so the surface gravity vanishes and the drawing uses $p = \arctan(u/m)$, $q = \arctan(v/m)$, shifted by $\pi$ from each region to the next, with the exteriors read from the isotropic chart and the interiors from Reissner-Nordström's own chart at $r_s = 2m$ and $r_q = m$.
Inside, $R_* = R + 2m\ln|R/m - 1| - m^2/(R - m) - m$ vanishes at $R = 0$, which puts the singularity on the vertical line $X = -\pi$.
The drawing's docstring lists the cells.
