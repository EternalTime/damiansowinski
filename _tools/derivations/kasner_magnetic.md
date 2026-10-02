# Kasner's universe with a magnetic field

The two charts of `kasner_magnetic` are written by `print_charts.py`, and this note records the source of each and what is checked.

## Sources

G. Rosen, "Symmetries of the Einstein-Maxwell equations", *Journal of Mathematical Physics* **3**, 313 (1962), and "Spatially homogeneous solutions to the Einstein-Maxwell equations", *Physical Review* **136**, B297 (1964), are the source of the solution; K. C. Jacobs, *The Astrophysical Journal* **155**, 379 (1969), calls it the pure-magnetic solution and credits Rosen with it in his section III.
D. Kastor and J. Traschen, "Melvin magnetic fluxtube/cosmology correspondence", *Classical and Quantum Gravity* **32**, 235027 (2015), arXiv:1507.05534, is the source of both charts as they are written here, and its equations are cited below from the arXiv version.
All four papers named by the candidate list, Rosen (1962), Thorne (1967) and Jacobs (1968, 1969), were verified on Crossref before anything was built, and Thorne's and Jacobs's second were read whole from the scans of the Astrophysics Data System.

## Kasner's time

Kastor and Traschen's generalized Melvin cosmology, the equation after their (23), is Harrison's transformation of Kasner's vacuum along $x^3$:

$$ds^2 = f^2\left[-dT^2 + (T/T_0)^{2p_1}(dx^1)^2 + (T/T_0)^{2p_2}(dx^2)^2\right] + f^{-2}(T/T_0)^{2p_3}(dx^3)^2, \qquad f = 1 + \frac{E^2}{4}(T/T_0)^{2p_3} .$$

The chart writes $b = E/2$ and reads the powers with $t$ in the unit $T_0$, as Kasner's own chart does.
They give the solution with an electric field along $z$ and note that duality trades it for the magnetic field $F_{xy}$ constant, with the metric unchanged.
In the chart $x^0 = ct$ the field is $F_{xy} = 2bp_3/c$ in units where it is an inverse length, so observers at rest measure $B = (2bp_3c/\sqrt{G})\,t^{p_3 - 1}/f^2$ in Gaussian units.
Only $p_3 > 0$ is listed: with $p_3 < 0$ the factor $f$ is $b^2t^{2p_3}(1 + t^{-2p_3}/b^2)$, the same family with its two ends exchanged.

On Kasner's circle the Ricci tensor is that of the field, with tension along $z$ and pressure across it,

$$R^t{}_t = R^z{}_z = -R^x{}_x = -R^y{}_y = -\frac{4b^2p_3^2\,t^{2p_3 - 2}}{c^2f^4}, \qquad R = 0,$$

and the Kretschmann scalar holds $p_3$ and $w = b^2t^{2p_3}$ alone, since $p_1p_2 = p_3^2 - p_3$ there:

$$K = \frac{16p_3^2\left[1 - p_3 + 6p_3(1 - p_3)w + 2(13p_3^2 - 2p_3 - 1)w^2 - 6p_3(3p_3 + 1)w^3 + (2p_3 + 1)(3p_3 + 1)w^4\right]}{c^4t^4f^8} .$$

At $b = 0$ it is Kasner's $16p_3^2(1 - p_3)/c^4t^4 = -16p_1p_2p_3/c^4t^4$.
At $p_3 = 1$ it is $64b^4(5 - 6b^2t^2 + 3b^4t^4)/c^4f^8$ with $t$ in its unit, finite at $t = 0$: the axisymmetric universe begins on a coordinate singularity, the edge of Milne's wedge, and for $p_3 < 1$ the scalar diverges as $t^{-4}$.
The Christoffel symbols and the Riemann tensor are printed for free exponents; the Ricci tensor, the scalars and the Einstein and Weyl tensors are printed as they stand on the circle through `OnSurface`, after `kasner_magnetic_check` holds each to the tensor of free exponents on the checker's parametrisation of the circle, the one Kasner's chart uses.
The stress of the field holds $t^{-2p_1 - 2p_2}$ where the curvature holds $t^{2p_3 - 2}$, so it is compared with the stated Ricci tensor after $p_1 = 1 - p_2 - p_3$, which makes the two one power.
`magnetic_powers` prints each value as one power of the time, with its whole exponent, times sums written in rising powers of $b^2t^{2p_3}$.

At late times $f \simeq b^2t^{2p_3}$, and in the proper time $\tilde t \propto t^{1 + 2p_3}$ the metric is Kasner's again with their (24),

$$\tilde p_1 = \frac{p_1 + 2p_3}{1 + 2p_3}, \qquad \tilde p_2 = \frac{p_2 + 2p_3}{1 + 2p_3}, \qquad \tilde p_3 = -\frac{p_3}{1 + 2p_3} .$$

## Rosen's chart

Kastor and Traschen's appendix writes Rosen's solutions as

$$ds^2 = -\frac{b^2\tan^{2(c_1 + c_2)}(t/2)}{\sin^4t}dt^2 + \frac{\tan^{2c_1}(t/2)}{\sin^2t}(dy^1)^2 + \frac{\tan^{2c_2}(t/2)}{\sin^2t}(dy^2)^2 + \sin^2t\,(dy^3)^2, \qquad c_1c_2 = 1,$$

with $p_1 = (c_1 - 1)p_3$, $p_2 = (c_2 - 1)p_3$ and $p_3 = 1/(c_1 + c_2 - 1)$.
The chart is the axisymmetric case $c_1 = c_2 = 1$, exponents $(0, 0, 1)$, where $\tan^2(t/2)/\sin^2t = 1/(1 + \cos t)^2$; it writes Rosen's time as $\eta$ and his length $b$ as $\ell$, since $t$ and $b$ are taken.
`kasner_magnetic_check` holds it to being Kasner's time at $(0, 0, 1)$ pulled back along $\tan(\eta/2) = bt$, with $x$ and $y$ halved, $z$ doubled and scaled by $b$, and $\ell = 2ct_0/b$.
Its field is $F_{xy} = 1/\ell$, so $B = (c^2/\sqrt{G})(1 + \cos\eta)^2/\ell$, finite at $\eta = 0$, Thorne's pancake, and $K = 8(1 + \cos\eta)^6(7\cos^2\eta + 2\cos\eta + 1)/\ell^4$.
The general case holds $\tan(\eta/2)$ to a power with a parameter in it, which the reader does not read, and is the first chart in another time.

## The drawings

The spacetime diagrams of Kasner's time and the embedding diagram are at $(p_1, p_2, p_3) = (-2/7, 3/7, 6/7)$, the point Kasner's own page draws, and $b = 1$, where the length along the field is greatest at $t = 1$ and the late exponents are $(10/19, 15/19, -6/19)$.
Across the field the cones are Kasner's, $x \pm (7/9)t^{9/7}$, and along it a ray keeps $z \pm (7t^{1/7} + (14/13)t^{13/7} + (7/25)t^{25/7})$.
In Rosen's chart a ray keeps $x \pm \ell\tan(\eta/2)$ or $z \pm (\ell/4)(\ln s + s^2 + s^4/4)$ with $s = \tan(\eta/2)$.
Rosen's views carry no moment of the embedding, which is drawn at other exponents.
There is no conformal diagram, as for Kasner's vacuum: a point of a plane of the time and one axis is a plane, and the axes differ.
