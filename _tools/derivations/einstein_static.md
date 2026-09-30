# The Einstein static universe

The three charts of `einstein_static.json` are written by `print_charts.py --metric einstein_static`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note derives what the charts are and what the diagrams draw.

## Step 1. The hyperspherical chart

The spatial slice is the round three sphere of radius $R$, the hypersurface $\xi_1^2 + \xi_2^2 + \xi_3^2 + \xi_4^2 = R^2$ of Euclidean space of four dimensions.
With $\xi_4 = R\cos\chi$ and $(\xi_1, \xi_2, \xi_3) = R\sin\chi\,(\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$ its metric is $R^2(d\chi^2 + \sin^2\chi\,d\Omega^2)$, so

$$ds^2 = -c^2dt^2 + R^2\left(d\chi^2 + \sin^2\chi\,d\theta^2 + \sin^2\chi\sin^2\theta\,d\phi^2\right),$$

with $\chi$ from the pole $\chi = 0$ to the antipode $\chi = \pi$.
Every component is taken in the chart $x^0 = ct$.

## Step 2. The curvature

The time direction is a flat factor, so every component of the Riemann tensor with a $t$ vanishes and the spatial ones are those of a space of constant curvature $1/R^2$:

$$R_{ijkl} = \frac{1}{R^2}\left(g_{ik}g_{jl} - g_{il}g_{jk}\right).$$

Contracting, $R_{ij} = 2g_{ij}/R^2$, $R_{tt} = 0$, the Ricci scalar is $6/R^2$ and the Kretschmann scalar is $12/R^4$, the same in every chart.
The Einstein tensor is $G^t{}_t = -3/R^2$ and $G^i{}_j = -\delta^i{}_j/R^2$.

The Weyl tensor vanishes: a product of a line with a space of constant curvature is conformally flat, which `print_charts.py` confirms component by component.

## Step 3. The field equations

With $G_{\mu\nu} + \Lambda g_{\mu\nu} = (8\pi G/c^4)T_{\mu\nu}$ and dust at rest, $T^t{}_t = -\rho c^2$ and $T^i{}_j = 0$, the two independent equations are

$$-\frac{3}{R^2} + \Lambda = -\frac{8\pi G\rho}{c^2}, \qquad -\frac{1}{R^2} + \Lambda = 0,$$

so $\Lambda = 1/R^2$ and $4\pi G\rho/c^2 = 1/R^2$.
This is Einstein's $\lambda = \kappa\rho/2 = 1/R^2$ with $\kappa = 8\pi G/c^2$, and the total mass is $M = 2\pi^2R^3\rho = \pi Rc^2/2G$.
No published component assumes these equations; the parameter's description states them.

## Step 4. The areal chart

With $r = R\sin\chi$ on the hemisphere $\chi < \pi/2$, $dr = R\cos\chi\,d\chi$ and $R^2d\chi^2 = R^2dr^2/(R^2 - r^2)$, so

$$ds^2 = -c^2dt^2 + \frac{R^2\,dr^2}{R^2 - r^2} + r^2d\Omega^2.$$

The chart ends at the equator $r = R$, where $g_{rr}$ diverges while the curvature stays $12/R^4$.
`print_charts.py` factors $R^2 - r^2$ into $(R + r)(R - r)$, and its `rewrite` puts it back whole, reading every rewritten value back against the one it replaces.

## Step 5. Einstein's projection

Einstein's own coordinates of 1917 are $\xi_1$, $\xi_2$ and $\xi_3$, the projection of the hemisphere $\xi_4 > 0$ onto its equatorial hyperplane, written $x$, $y$ and $z$.
Eliminating $\xi_4 = \sqrt{R^2 - x^2 - y^2 - z^2}$,

$$d\xi_4^2 = \frac{(x\,dx + y\,dy + z\,dz)^2}{R^2 - x^2 - y^2 - z^2},$$

which is his $\gamma_{\mu\nu} = \delta_{\mu\nu} + x_\mu x_\nu/(R^2 - \rho^2)$.
Its Christoffel symbols are $\Gamma^i{}_{jk} = x^i\gamma_{jk}/R^2$, so every geodesic equation reads

$$\ddot{x}^i + \frac{x^i}{R^2}\left(\dot{x}^2 + \dot{y}^2 + \dot{z}^2 + \frac{(x\dot{x} + y\dot{y} + z\dot{z})^2}{R^2 - x^2 - y^2 - z^2}\right) = 0,$$

which `print_charts.py` writes by hand for this chart and for the areal chart and checks against the equations it prints from the Christoffel symbols.

## Step 6. The spacetime diagrams

On the plane of $t$ and $\chi$ the metric is $-c^2dt^2 + R^2d\chi^2$, so the rays conserve $ct \pm R\chi$, which `null_rays.py --verify` measures.
Along the great circle through the pole the view through the pole draws $\phi = 0$ on the right and $\phi = \pi$ on the left, and its two edges $x = \pm\pi$ are the one world line of the antipode.
In the areal chart $dt/dr = \pm R/(c\sqrt{R^2 - r^2})$, so the rays conserve $ct \mp R\arcsin(r/R)$ and reach $r = R$ after $\pi R/2c$.

## Step 7. The conformal diagram

With $\eta = ct/R$ the plane of $t$ and $\chi$ is $R^2(-d\eta^2 + d\chi^2)$, already conformally flat, and $p, q = (\eta \mp \chi)/2$ draw it as the strip $0 \le \chi \le \pi$ with $X = \chi$ and $T = \eta$.
The other three maps are the ones the collection's own conformal diagrams use, and each has $X = \chi$ and $T = \eta$ of this strip:

- Minkowski space, $ct \pm r = R\tan((\eta \pm \chi)/2)$, with metric $\Omega^{-2}$ times this one, $\Omega = 2\cos((\eta + \chi)/2)\cos((\eta - \chi)/2)$, on $|\eta| + \chi < \pi$;
- de Sitter space of radius $\ell = R$, with metric $1/\cos^2\eta$ times this one, on $|\eta| < \pi/2$;
- anti-de Sitter space of radius $L = R$, $r = L\tan\chi$, with metric $1/\cos^2\chi$ times this one, on $\chi < \pi/2$.

`conformal.py` checks at random points of each that the published metric of the Einstein static universe, pulled back through the map, is $\Omega^2$ times the other spacetime's published metric on its plane of $t$ and $r$, with the $\Omega^2 = \sin^2\chi/g_{\theta\theta}$ that the sphere through the point requires.

## Step 8. The embedding diagram

The equator of one moment of the hyperspherical chart is $R^2(d\chi^2 + \sin^2\chi\,d\phi^2)$: $\rho = R\sin\chi$, $d\rho/d\chi = R\cos\chi$ and $dz/d\chi = R\sin\chi$, the sphere $z = -R\cos\chi$ of radius $R$, the same at every moment.
The areal chart's slice, $R^2dr^2/(R^2 - r^2) + r^2d\phi^2$, is checked to be its lower hemisphere, $z = -\sqrt{R^2 - r^2}$.

## Step 9. The regions on the sphere

At the conformal time $\eta$ each spacetime drawn inside the strip covers the values of $\chi$ its map reaches from its own chart, and `conformal.covered()` finds them by bisection in that chart's time at every radius.
Minkowski space, through $ct \pm r = R\tan((\eta \pm \chi)/2)$, covers $\chi < \pi - |\eta|$, all of the sphere but the antipode, spatial infinity $i^0$, at $\eta = 0$.
De Sitter space enters by its closed slicing, $X_0 = \ell\sinh(ct/\ell)$ on its hyperboloid, whose global coordinates put $\tan\eta = X_0/\ell$ with $\chi$ unchanged, so it covers the whole sphere for $|\eta| < \pi/2$ and none of it beyond; its static patch alone covers only $|\eta| + \chi < \pi/2$.
Anti-de Sitter space, through $r = L\tan\chi$ and $ct = L\eta$, covers the hemisphere $\chi < \pi/2$ at every $\eta$.
These are Hawking and Ellis's regions, and `embedding.py` checks them at $\eta$ from $-3$ to $3$ in steps of $1/2$, while `conformal.py` checks the region each view draws in the strip against the same map.
The embedding diagram's sphere is the moment $t = 0$, $\eta = 0$, the slice both diagrams draw, and it is shaded there while the conformal diagram shows the region.

## Checking

    <venv>/bin/python _tools/derivations/print_charts.py --metric einstein_static
    <venv>/bin/python _tools/derivations/verify_metrics.py --system einstein_static/hyperspherical --system einstein_static/static_areal --system einstein_static/einstein_cartesian
    <venv>/bin/python _tools/derivations/embedding.py --metric einstein_static
    <venv>/bin/python _tools/derivations/null_rays.py --metric einstein_static
    <venv>/bin/python _tools/derivations/conformal.py --metric einstein_static
    python3 _tools/build_mfs_data.py
