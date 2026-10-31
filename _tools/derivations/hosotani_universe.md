# Hosotani universe

## Sources and conventions

Yutaka Hosotani, *Exact solution to the Einstein-Yang-Mills equation*, Physics Letters B 147 (1984), 44-46, DOI 10.1016/0370-2693(84)90588-4, establishes the closed Friedmann solution and its gauge instability in the publisher-supplied abstract.
The original full text remains unread: the publisher route refused access, and its advertised text-mining route returned HTTP 400.
G. W. Gibbons and Alan R. Steif, *Yang-Mills cosmologies and collapsing gravitational sphalerons*, Physics Letters B 320 (1994), 245-252, DOI 10.1016/0370-2693(94)90652-1, arXiv:hep-th/9311098, supply readable primary equations.
Their section 2 explicitly identifies the stationary magnetic configuration with Hosotani through reference 7.
We use signature $(-,+,+,+)$ and units $c=\hbar=1$.
The gauge coupling $e$ is dimensionless and $G$ has dimensions of length squared in these units.

## Homogeneous gauge field

Let $U=x^4\mathbf{1}+ix^a\tau_a$ map the unit spatial 3-sphere to $SU(2)$, where $\sum_{A=1}^4(x^A)^2=1$ and $\tau_a$ are Pauli matrices.
Define the left invariant forms $\omega^a=-i\operatorname{Tr}(\tau_aU^{-1}dU)$, with $d\omega^a=\tfrac12\epsilon_{abc}\omega^b\wedge\omega^c$ and $d\Omega_3^2=\tfrac14\sum_a(\omega^a)^2$.
These are Gibbons and Steif's equations (2.2), (2.3), and (2.6).
Their connection and curvature, equations (2.7) and (2.8), are

$$A=\frac{i}{e}f(\eta)U^{-1}dU=-\frac{f}{e}\omega^a\frac{\tau_a}{2},\qquad F=dA-ieA\wedge A=F^a\frac{\tau_a}{2},$$
$$F^a=-\frac{f'}{e}d\eta\wedge\omega^a+\frac{f(f-1)}{2e}\epsilon_{abc}\omega^b\wedge\omega^c.$$

The commutator in $F$ is essential to the magnetic term and to the stationary nonzero solution.
A right translation rotates the forms in the adjoint representation and is compensated by a global gauge rotation; a left translation leaves them invariant.
The stress tensor is invariant under the full spatial rotation group.
This is one $SU(2)$ connection, with its internal directions coupled by the Lie bracket.

Equation (2.10) reduces the Yang-Mills equation to

$$f''=-\frac{dV}{df},\qquad V=2(f^2-f)^2,\qquad E=\frac12(f')^2+V.$$

Hosotani's stationary gauge configuration is $f=\tfrac12$, with $E=\tfrac18$.
The stationary solutions $f=0$ and $f=1$ have zero curvature and zero energy.
The time dependent solutions at other $E$ form a broader homogeneous gauge family and are not assumed as additional Hosotani branches.
For $f=\tfrac12+\delta f$, the homogeneous linearized equation is $\delta f''=2\delta f$.
This verifies a growing homogeneous gauge perturbation, while the count of all unstable gauge harmonics rests on Hosotani's abstract until the original analysis is read.

## Stress tensor and scale factor

Gibbons and Steif's equation (2.12) uses

$$T_{\mu\nu}=\sum_a\left(F^a_{\mu\alpha}F^a_\nu{}^\alpha-\frac14g_{\mu\nu}F^a_{\alpha\beta}F^{a\alpha\beta}\right).$$

For the orthonormal spatial coframe $\vartheta^a=(a/2)\omega^a$, the electric and magnetic components are $E^a_i=-2f'\delta^a_i/(ea^2)$ and $B^a_i=4f(f-1)\delta^a_i/(ea^2)$.
Thus $\rho=12E/(e^2a^4)$ and $p=\rho/3$, agreeing with the paragraph following (2.12).
At $f=\tfrac12$, $E^a_i=0$, $B^a_i=-\delta^a_i/(ea^2)$, and $\rho=3/(2e^2a^4)$.
The word stationary refers to the connection in conformal coordinates; the magnetic strength measured in an orthonormal frame varies as $a^{-2}$.

The vanishing cosmological constant and the closed radiation Friedmann equation imply

$$\dot a^2+1=\frac{\ell^2}{a^2},\qquad (a')^2+a^2=\ell^2,\qquad \ell^2=\frac{4\pi G}{e^2}.$$

Here a dot differentiates cosmic time $t$ and a prime differentiates dimensionless conformal time $\eta$.
The primary solution after (2.12), $a=\sqrt{32\pi E}\,\sqrt{G}\,e^{-1}\sin\eta$, specializes to

$$a(\eta)=\ell\sin\eta,\qquad t=\ell(1-\cos\eta),\qquad a(t)=\sqrt{t(2\ell-t)}.$$

Expansion occupies $0<\eta<\pi/2$ and contraction occupies $\pi/2<\eta<\pi$.
The maximum radius is $\ell$ at $t=\ell$; the total cosmic lifetime is $2\ell$.
These two phases belong to one connected solution, with a regular turnaround.
The endpoints at $t=0$ and $t=2\ell$ are singular.
No cosmological constant extension or static gravitational branch is inferred from a generic FLRW equation.

## Coordinate charts

Gibbons and Steif's equation (2.1) supplies cosmic time and conformal time charts with the round spatial metric on $S^3$.
Using $d\Omega_3^2=d\chi^2+\sin^2\chi(d\theta^2+\sin^2\theta\,d\phi^2)$ makes their charts explicit:

$$ds^2=-dt^2+t(2\ell-t)d\Omega_3^2,$$
$$ds^2=\ell^2\sin^2\eta(-d\eta^2+d\Omega_3^2).$$

We take $0<t<2\ell$, $0<\eta<\pi$, $0<\chi<\pi$, $0<\theta<\pi$, and $0\leq\phi<2\pi$ on each angular coordinate patch.
The coordinate degeneracies at the polar axes are those of the round 3-sphere.
The areal coordinate $r=\sin\chi$ covers a hemisphere and is not used as a global spatial coordinate.

## Drawings

The spacetime views use cosmic and conformal time with the polar distance $\chi$.
The radial null equation is $d\chi/d\eta=\pm1$.
The conformal orbit space is the rectangle $0<\eta<\pi$, $0\leq\chi\leq\pi$, with spacelike singularities at its past and future edges and regular symmetry centres at its sides.
At fixed cosmic time the equatorial spatial 2-sphere has the Euclidean embedding $X^2+Y^2+Z^2=a(t)^2$.
A single view at maximum expansion, $t=\ell$, suffices to represent that intrinsic geometry.
The curvature figure reads only the stored Kretschmann scalar over $0.02\ell\leq t\leq1.98\ell$, on a logarithmic scale from $10$ to $10^8$ in $K\ell^4$.
The Ricci scalar vanishes throughout the spacetime.
Each drawing uses $\ell=1$ as its unit; the gauge source normalization remains fixed.

## Checks

Both source-specialized charts are generated by `hosotani_universe.chart_spec`, registered in `print_charts.py`, with point comparisons of the Einstein source, the chart pullback, and both curvature scalars.
`test_hosotani_universe.py` checks the gauge curvature directly with Pauli matrices and the Maurer-Cartan equation, the isotropic stress, the source normalization $\ell^2=4\pi G/e^2$, both Friedmann equations, the growing homogeneous gauge mode, and the coordinate transformation of the stored metric components.
These source checks remain distinct from the general tensor verification.
