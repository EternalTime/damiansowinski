# The Hawking-Moss instanton

## Sources

Stephen Hawking and Ian G. Moss, *Supercooled phase transitions in the very early universe*, Physics Letters B 110 (1982), 35-38, DOI 10.1016/0370-2693(82)90946-7.
The founding journal metadata is verified; its original full text remains unread.
Ian G. Moss's primary self-author line in arXiv:2303.01148, p. 1, resolves the inconsistent I.L. initials in the founding metadata.
No quotation or equation is taken from the unread original.

Erick J. Weinberg, *Hawking-Moss Bounces and Vacuum Decay Rates*, Physical Review Letters 98 (2007), 251303, arXiv:hep-th/0612146v2, pp. 1-4, supplies the radius and exponent in equation (1), the thermal expansion in equation (2), and the redshifted energy in footnote 5.
Adam R. Brown and Erick J. Weinberg, *Thermal derivation of the Coleman-De Luccia tunneling prescription*, Physical Review D 76 (2007), 064003, arXiv:0706.1573v3, pp. 1-28, supplies the coordinate systems in equations (3.6)-(3.7), (5.3)-(5.7), and footnote 7.
Both later primary texts were read in full.
Daiki Saito and Naritaka Oshita, JHEP 02 (2025), 187, arXiv:2409.03978v3, was also read in full; its remote-transition probability assumes a homogeneous scalar on fixed de Sitter space, subdominant kinetic energy, small potential variation, and no thermal or stochastic noise.
Its linear-potential result supplies no general decay future.

## The homogeneous saddle

We take a canonical scalar with a positive false vacuum energy density $V_F$ and a higher positive local maximum $V_{\mathrm{top}}$.
The saddle holds $\phi = \phi_{\mathrm{top}}$ with $V'(\phi_{\mathrm{top}}) = 0$ and $V''(\phi_{\mathrm{top}}) < 0$.
In units $c = \hbar = 1$, its radius obeys $\ell^2 = 3/(8\pi G V_{\mathrm{top}})$.
The false vacuum sphere has its own radius $\ell_F^2 = 3/(8\pi G V_F)$.

## Euclidean coordinate charts

The proper-distance form is $ds^2 = d\xi^2 + \ell^2\sin^2(\xi/\ell)d\Omega_3^2$, for the whole interval $0 \le \xi \le \pi\ell$.
Weinberg p. 2 states the sphere; Brown and Weinberg equation (5.6) uses a dimensionless radial angle and a shifted 3-sphere latitude.
The rescalings are $\xi_{\mathrm{there}} = \xi/\ell$ and $\tau'_{\mathrm{there}} = \chi - \pi/2$.

Brown and Weinberg equations (3.6)-(3.7) write the same sphere as $ds^2 = f\,d\tau^2 + dr^2/f + r^2d\Omega_2^2$, $f = 1-r^2/\ell^2$.
The Euclidean time has period $2\pi\ell$; near $r=\ell$, a proper polar radius $\varrho$ gives $d\varrho^2 + \varrho^2d(\tau/\ell)^2$, fixing that period by smoothness.
Equation (5.3) replaces $r$ by $\ell\sin\eta$, giving $ds^2 = \cos^2\eta\,d\tau^2 + \ell^2(d\eta^2 + \sin^2\eta\,d\Omega_2^2)$.
At $r=0$ the 2-sphere contracts; at $r=\ell$ the Euclidean time circle contracts.
The endpoints are regular polar degeneracies of these coordinates.

## Lorentzian continuation

Brown and Weinberg equations (1.3) and (5.4) continue $\tau = it$ to the static patch, $ds^2 = -f\,dt^2 + dr^2/f + r^2d\Omega_2^2$.
Equation (5.7) continues $\chi = \pi/2 + i\psi$ in the sphere chart to $ds^2 = d\xi^2 + \ell^2\sin^2(\xi/\ell)(-d\psi^2 + \cosh^2\psi\,d\Omega_2^2)$.
Their footnote 7 continues at the sphere's equator, $\xi = \pi\ell/2 + it$, giving $ds^2 = -dt^2 + \ell^2\cosh^2(t/\ell)d\Omega_3^2$.
All three are exact constant-field de Sitter geometries.
The field remains at the unstable maximum until perturbed; its rolling evolution requires a potential and initial fluctuations beyond these exact saddles.

The coordinate checks pull every metric back from a round sphere in Euclidean five-space or a de Sitter hyperboloid in Minkowski five-space at three regular points and three different radii.
They also compare $R_{ab} = 3g_{ab}/\ell^2$ at those points.
Every tensor is generated from the corresponding line element and read back by the chart printer.

## Action and thermal activation

The Euclidean action is

$$S_E = \int d^4x\sqrt{g}\left[-\frac{R}{16\pi G} + \frac12(\nabla\phi)^2 + V(\phi)\right].$$

For the constant scalar, $R_{ab} = 8\pi GVg_{ab}$ and $R = 32\pi GV = 12/\ell^2$.
The on-shell integrand is $-V$, and the volume of a 4-sphere is $8\pi^2\ell^4/3$.
Consequently $S_E = -\pi\ell^2/G = -3/(8G^2V)$.
Subtracting the false vacuum action gives

$$B = \frac{3}{8G^2}\left(\frac{1}{V_F}-\frac{1}{V_{\mathrm{top}}}\right),\qquad e^{-B}.$$

This agrees with Weinberg equation (1); it specifies only the leading exponential weight.
Writing $\Delta V = V_{\mathrm{top}}-V_F$ and expanding for $\Delta V/V_F \ll 1$ gives $B \simeq 8\pi^2\ell_F^4\Delta V/3$.
The static observer at the centre has temperature $T_F = 1/(2\pi\ell_F)$.
Its Killing energy difference is

$$\Delta E = \int_0^{\ell_F}4\pi r^2\,dr\,\frac{\sqrt{1-r^2/\ell_F^2}}{\sqrt{1-r^2/\ell_F^2}}\,\Delta V
= \frac{4\pi}{3}\ell_F^3\Delta V,$$

so $\Delta E/T_F$ agrees with the leading expansion of $B$.
The lapse in the numerator is essential: multiplying $\Delta V$ by the unweighted proper static volume gives a different energy.
Weinberg footnote 5 and Brown and Weinberg equation (3.4) supply this energy integral.
The fixed geometry thermal derivation concerns one horizon volume and leaves the state beyond its horizon unrestricted.

The scalar quadratic operator is $-\nabla_{S^4}^2 + V''(\phi_{\mathrm{top}})$.
Its sphere harmonic eigenvalues are $n(n+3)/\ell^2 + V''(\phi_{\mathrm{top}})$ for integers $n\ge0$.
A local crossing of one eigenvalue supplies no global criterion for existence or dominance of all Coleman-De Luccia saddles.
Brown and Weinberg footnote 6 and the later primary counterexamples recorded in the source ledger establish this qualification.
The O(4)-invariant sector has its own multiplicities and must be distinguished from the full sphere harmonic degeneracies.

## Curvature and embeddings

The generated invariants are $R = 12/\ell^2$ and $K = 24/\ell^4$.
The curvature figure reads these stored expressions through curvature.py; it does not derive curvature again.
At $\chi = \theta = \pi/2$, the Euclidean meridian metric is $d\xi^2 + \ell^2\sin^2(\xi/\ell)d\phi^2$.
The isometric surface of revolution has $\rho = \ell\sin(\xi/\ell)$ and $z = \ell\cos(\xi/\ell)$.
The full interval $0\le\xi\le\pi\ell$ retains both hemispheres and both regular poles.
All drawings use $\ell=1$.

## Diagram domains

NOT_DRAWN, four_sphere, thermal, thermal_angle: these metrics are positive definite at regular points, so the only null tangent is zero and no light cone or causal conformal diagram exists.
The three Lorentzian light-cone views depict only the specified constant-field mathematical continuations.
The conformal view is restricted to the static continuation's patch, with its regular centre and cosmological horizons.
It describes no vacuum-decay wall, rolling field, or global decay future.
The Euclidean meridian surface has no Lorentzian moment to mark on those continuation diagrams; their slice exception records that distinction.
All figure windows remain within aspect ratios 1:2 and 2:1, with margins for the entire geometry through rotation.
