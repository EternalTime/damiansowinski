# The Klebanov-Witten throat and its conifold

## Sources and coordinates

Igor R. Klebanov and Edward Witten's hep-th/9807080v3, equations (1), (3), (4), and (10), fixes the cone, its Einstein base, and the throat.
Their equation (4) in the near-horizon limit has

$$ds_{10}^2=\frac{r^2}{L^2}(-c^2dt^2+dx^2+dy^2+dz^2)+\frac{L^2}{r^2}dr^2+L^2ds_{T^{1,1}}^2.$$

We use the angular names $(\theta,\phi,\alpha,\beta,\psi)=(\theta_1,\phi_1,\theta_2,\phi_2,\psi)$ of their equation (10).
The first four angles describe two 2-spheres, and $\psi$ has period $4\pi$.
The time is unwrapped; $r>0$ covers the Poincaré patch.
The singular angular poles are the familiar degeneracy of spherical coordinates.

The six-dimensional transverse space is a Riemannian cone,

$$ds_6^2=dr^2+r^2ds_{T^{1,1}}^2.$$

Klebanov and Witten attribute its Ricci flat Kähler metric to Philip Candelas and Xenia C. de la Ossa's 1990 work.
The publisher refused access to that full text; Crossref verifies its primary bibliographic record, and the explicit equation and attribution are open in Klebanov and Witten, pp. 4-5.
Larry James Romans's own Caltech thesis, Chapter 2, section 2.2, reproduces his 1985 compactification on the $p=q=1$ Einstein coset.
Its provenance footnote on p. 10 and reference 1 on p. 29 connect the derivation to his Physics Letters article.
The journal article's full text remains an access debt.

The hypersurface $\sum_{i=1}^4 z_i^2=0$, its Kähler potential, and its quotient construction in equations (7)-(15) describe the transverse cone through constrained complex variables.
They do not constitute independent real charts of the ten-dimensional throat.
The conformally related cylinder of equation (2) belongs to the cone construction.
The resolved and deformed conifolds are different Riemannian metrics.
The Poincaré chart is the explicit sourced spacetime coordinate system; its complete derived tensors are retained.

## Curvature

For the normalized base $h$, direct angular differentiation gives $\operatorname{Ric}(h)=4h$, $R(h)=20$, and $|\operatorname{Riem}(h)|^2=136$.
For its cone, $R_{rr}=R_{ra}=0$ and $R_{ab}=\operatorname{Ric}(h)_{ab}-4h_{ab}=0$.
Its squared Riemann tensor is $(136-2\cdot20)/r^4=96/r^4$; the cone is singular at the tip.
The five-dimensional anti-de Sitter factor has Ricci eigenvalue $-4/L^2$, scalar $-20/L^2$, and squared Riemann tensor $40/L^4$.
The compact factor has Ricci eigenvalue $4/L^2$, scalar $20/L^2$, and squared Riemann tensor $136/L^4$.
The throat's scalar is zero, its Ricci tensor is nonzero, and its squared Riemann tensor is $(40+136)/L^4=176/L^4$.
The direct product has no mixed Riemann components between its two factors.
All these statements follow from coordinate differentiation, with the base's squared Riemann tensor also evaluated at three independent pairs of polar angles.

## Drawings and their scope

The spacetime diagram is the induced $(ct,r)$ plane at $x=y=z=0$, $\theta=\alpha=\pi/2$, and $\phi=\beta=\psi=0$.
Its radial null geodesics satisfy $dr/d(ct)=\pm r^2/L^2$; their constants are $ct\pm L^2/r$.
This plane fixes the internal point and the three brane coordinates, and its curves are null geodesics of the full product.
It does not encode the trajectories carrying internal angular momentum or displacement along the branes.

The conformal diagram compactifies the same two-dimensional induced plane through the Poincaré map of anti-de Sitter space of two dimensions.
It shows that plane's timelike boundary and Poincaré horizon, with its slice named in the caption.
It does not describe the full ten-dimensional causal structure or identify the other eight dimensions with single events.

The embedding is the spatial surface $(r,\phi)$ at $t=x=y=z=0$, $\theta=\alpha=\pi/2$, and $\beta=\psi=0$.
Its metric is $L^2dr^2/r^2+(L^2/6)d\phi^2$, so it embeds as a cylinder of radius $L/\sqrt6$ and height $L\ln(r/L)$.
The cylinder extends without end in both directions; its finite reach for drawing is $L/2\leq r\leq5L/2$.
It resolves one circle of the compact base and the throat's radial direction, not the full spatial geometry of nine dimensions.
The moment $t=0$ on this cylinder is registered on the radial null and conformal diagrams with view id throat, and those slices have the same radial endpoints.

A complete isometric embedding of nine spatial dimensions in Euclidean three-space is NOT_DRAWN because the ambient dimension is too small.
A separate embedding of the six-dimensional transverse cone is NOT_DRAWN because it is a different Riemannian space, not a moment of the throat; the cone's metric and singular curvature are derived explicitly above.
The required throat embedding is the induced cylindrical surface.
No chart or required kind of drawing is omitted.

The curvature drawing must display the cancellation between the two Ricci sectors and the zero scalar, from the stored curvature.
A radial profile cannot distinguish homogeneous sectors, since neither scalar changes with position.

## Verification

The chart check compares the sourced product and its opposite Ricci eigenvalues at three regular points, without a symbolic determinant or symbolic simplification inside the check.
The dedicated physical checks also compute the base and cone directly, read the throat's derived tensors, and compare its scalar invariants and the curvature decomposition.
The spacetime and conformal slices carry the embedding's moment version and view id.
