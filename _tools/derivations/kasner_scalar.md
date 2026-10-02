# Kasner's universe with a scalar field

The three charts of `kasner_scalar` are written by `print_charts.py`, and this note records the source of each and what is checked.

## Sources

V. A. Belinskii and I. M. Khalatnikov, "Effect of scalar and vector fields on the nature of the cosmological singularity", *Soviet Physics JETP* **36**, 591 (1973), submitted 27 April 1972, is the source of the solution, and its equations are cited below by their numbers.
K. C. Jacobs, "Spatially homogeneous and Euclidean cosmological models with shear", *The Astrophysical Journal* **153**, 661 (1968), has the same metric as the universe of a fluid with pressure equal to energy density, his (49).
T. Damour, M. Henneaux and H. Nicolai, "Cosmological billiards", *Classical and Quantum Gravity* **20**, R145 (2003), is the source of the logarithmic time, their (3.20) to (3.28).

## The synchronous chart

Their (2.6) to (2.8), in units where $c$ and the Einstein gravitational constant are one:

$$ds^2 = -dt^2 + t^{2p_1}dx^2 + t^{2p_2}dy^2 + t^{2p_3}dz^2, \qquad \varphi = q\ln t, \qquad \sum p_i = 1, \quad \sum p_i^2 = 1 - q^2 .$$

Their field equations are $R_{ik} = \varphi_{;i}\varphi_{;k}$, (2.5), so with $c$ and $G$ kept the field is $\varphi = (qc^2/\sqrt{8\pi G})\ln t$ and $R_{\mu\nu} = (8\pi G/c^4)\,\partial_\mu\varphi\,\partial_\nu\varphi$.
Their (2.9) bounds $q^2 \le 2/3$, and all three exponents are positive wherever $q^2 > 1/2$: an exponent that is not positive leaves the other two summing to at least one, so their squares sum to at least $1/2$.
Jacobs writes the exponents as $p_i = 1/3 + (2\delta/3)\sin(\psi + 2\pi(i-1)/3)$ with $0 \le |\delta| < 1$, which is the same surface with $q^2 = (2/3)(1 - \delta^2)$, and his density, $(1 - \delta^2)/24\pi t^2$ in units $G = c = 1$, is the $q^2c^2/16\pi Gt^2$ of the History.

For free exponents the Ricci tensor is

$$R_{tt} = \frac{\sum p_i - \sum p_i^2}{c^2t^2}, \qquad R_{ii} = \frac{p_i\left(\sum p_j - 1\right)t^{2p_i - 2}}{c^2},$$

so the two conditions are the field equations, as they are for Kasner's vacuum in `kasner.md`.
The Christoffel symbols and the Riemann tensor are printed for free exponents, which hold on the surface too.
The Ricci tensor, both scalars and the Einstein and Weyl tensors are printed as they stand on the surface, written in the exponents and $q$: `OnSurface` in `print_charts.py` hands them to the printer, and `kasner_scalar_check` holds each to the tensor of free exponents on the parametrisation below before anything is written.
On the surface

$$R = -\frac{q^2}{c^2t^2}, \qquad K = \frac{3q^4 - 16p_1p_2p_3}{c^4t^4}, \qquad C_{txtx} = \frac{\left(3p_1(1 - p_1) - q^2\right)t^{2p_1 - 2}}{3c^2}, \qquad C_{xyxy} = \frac{\left(6p_1p_2 - q^2\right)t^{2p_1 + 2p_2 - 2}}{6c^2},$$

and the Weyl tensor vanishes at the isotropic point $p_i = 1/3$, $q^2 = 2/3$, the flat Friedmann universe of a stiff fluid.
The Kretschmann scalar follows from $K = 4\left[\sum p_i^2(p_i - 1)^2 + \sum_{i<j}p_i^2p_j^2\right]/c^4t^4$ and the elementary symmetric functions $e_1 = 1$, $e_2 = q^2/2$ and $e_3 = p_1p_2p_3$.

## The surface, as the checker holds it

`PARAMETER_RELATIONS` in `verify_metrics.py` needs a rational parametrisation.
The surface is a sphere about the isotropic point in the plane $\sum p_i = 1$, with Kasner's circle for its equator $q = 0$.
A point of it is the point $u$ of Kasner's circle drawn toward the isotropic point by the factor $r = (3v^2 - 2)/(3v^2 + 2)$:

$$p_i = \tfrac{1}{3} + r\left(k_i(u) - \tfrac{1}{3}\right), \qquad q = \frac{4v}{3v^2 + 2},$$

with $k_i(u)$ the parametrisation of `kasner.md`.
Then $\sum p_i^2 = 1/3 + 2r^2/3$ and $q^2 = (2/3)(1 - r^2)$.

## The logarithmic chart

With $\tau = -\ln(t/t_0)$ and $\ell = ct_0$,

$$ds^2 = -\ell^2e^{-2\tau}d\tau^2 + e^{-2p_1\tau}dx^2 + e^{-2p_2\tau}dy^2 + e^{-2p_3\tau}dz^2, \qquad \varphi \propto -q\tau .$$

This is the gauge of Damour, Henneaux and Nicolai, whose lapse is the root of the determinant of the spatial metric: their (3.21), $\beta^\mu = v^\mu\tau + \beta_0^\mu$, with the metric $-N^2d\tau^2 + \sum e^{-2\beta^i}dx_i^2$ and $\sum v^i = 1$.
The chart is checked to be the synchronous one pulled back along $t = t_0e^{-\tau}$, the powers read with $t$ in the unit $t_0$.

## The chart of five dimensions

Their (4.5) and (4.6): Kasner's vacuum of five dimensions,

$$ds^2 = -c^2dT^2 + T^{2s_1}dx^2 + T^{2s_2}dy^2 + T^{2s_3}dz^2 + T^{2s_5}dw^2, \qquad \sum s_a = \sum s_a^2 = 1 .$$

By their (3.21), (3.28) and fourth footnote the scalar field is the size of the fifth dimension and the metric of four dimensions is $\sqrt{g_{ww}}$ times the first four terms.
Its proper time is $t = T^{1 + s_5/2}/(1 + s_5/2)$, so

$$p_i = \frac{2s_i + s_5}{2 + s_5}, \qquad q = \frac{\sqrt{6}\,s_5}{2 + s_5}, \qquad s_5 = \frac{2q}{\sqrt{6} - q}, \qquad s_i = \frac{\sqrt{6}\,p_i - q}{\sqrt{6} - q},$$

the last two as their footnote prints them.
`kasner_scalar_check` holds $\sum p_i = 1$ and $\sum p_i^2 = 1 - q^2$ on the surface of five dimensions, whose parametrisation is the second point where the line from $(1, 0, 0, 0)$ along $(-(a + b + 1), a, b, 1)$ meets the sphere.
`REGIONS` in `metric_tags.py` counts the two charts of four dimensions, since the spacetime is the one with the scalar field in it.

## The drawings

Every drawing is at $(p_1, p_2, p_3) = (2/13, 4/13, 7/13)$ and $q = 10/13$, a rational point of the surface with all three exponents positive, and the chart of five dimensions at the exponents that reduce to it, about $(-0.234, -0.009, 0.327, 0.916)$.
There is no conformal diagram, as for Kasner's vacuum: a point of a plane of the time and one axis is a plane, and the three axes differ.
