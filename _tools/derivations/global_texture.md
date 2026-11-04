# Collapsing global texture

## Sources and approximation

Nötzold's Fermilab preprint FERMILAB-Pub-90/64-A, equation (10), printed p. 4, supplies the first-order proper-time and proper-distance metric.
Durrer, Heusler, Jetzer, and Straumann's author preprint ZU-TH 9/91, equations (0.30) and (0.46), printed pp. 6 and 8, supplies the independently matching angular factor.
The published journal citations are in the metric file; equation numbers here refer to the fetched preprints.
The two History readings approved the page, its fields, relations, references, nonlinear equations, and drawing plan through steering messages 002 and 003.

## Verified chart

The collapse chart has $t<0$, $r>0$, and $ds^2=-dt^2+dr^2+r^2H(t,r)d\Omega^2$.
Here $H=1-2\epsilon+2\epsilon(t/r)\arctan(r/t)$ and $\epsilon=8\pi G\eta^2$ in units $c=\hbar=1$.
This line element is retained only through first order as a texture solution.
The published tensors are exact tensors of the truncated metric, including the higher powers introduced by its inverse, and are not claimed to solve the full Einstein-sigma equations at finite epsilon.

## Numerical gauges and continuation

The nonlinear areal-radius AB gauge in equations (0.6) and (0.18)-(0.23) is not published as a chart here: its numerical metric functions and shooting solution have not been independently reproduced and verified.
The nonlinear omega gauge in equations (0.30)-(0.38) is stated as a system of equations on the page, without substituting first-order tensors for its full solution.
The coordinate transformation between those full gauges, equations (0.42)-(0.45), depends on their numerical functions.
No chart across the collapse is published.
Durrer and colleagues' section 3.2 distinguishes the even first-order omega gauge from the extension regular near $t=0$, and their section 3.3 states the sigma model's small-core limitation.
Nötzold's equations (7) describe multiple post-collapse continuations, so no continuation is selected here.

## Drawings

All four drawings describe the pre-collapse chart at $\epsilon=1/20$, with an arbitrary length $\ell$ as unit.
The radial null rays satisfy $dt/dr=\pm1$ exactly because the radial metric is $-dt^2+dr^2$.
Their plane covers $0\le r/\ell\le4$ and $-4\le t/\ell\le-0.05$, and marks the self-similar incoming light cone $r=-t$.

The conformal coordinates are $p=\arctan((t-r)/\ell)$ and $q=\arctan((t+r)/\ell)$, with displayed $X=q-p$ and $T=q+p$.
The radial metric is $-\ell^2\sec^2p\sec^2q\,dp\,dq$, a positive conformal factor times a flat null plane.
The pre-collapse region is the triangle with vertices $(X,T)=(0,-\pi),(0,0),(\pi,0)$.
Its future boundary is $t=0$, and the diagonal boundary from past timelike infinity to spatial infinity is past null infinity.
No future continuation is drawn, since the source model requires a choice of post-collapse continuation and does not resolve the core event.
The conformal program independently checks the coordinate map against the stored radial metric.

The embedding is the equator at $t=-\ell$, from $r=0$ to $r=4\ell$.
Its Euclidean cylindrical radius is $\rho=\sqrt{g_{\phi\phi}}=r\sqrt H$, and $dz/dr=\sqrt{1-(d\rho/dr)^2}$.
The shared embedding program takes that derivative from the stored metric and integrates the height by adaptive quadrature.
Its circumference, radial lengths and interpolated profile pass the shared isometry checks.
At the centre $H=1-2\epsilon r^2/(3t^2)+O(r^4/t^4)$, so the surface closes regularly and $d\rho/dr\to1$.
The surface has time $-1$ in units of $\ell$, label $t=-\ell$, and view id `collapse`.
Both causal diagrams register this same view, surface 0, from the same radial endpoints 0 and 4.
The null-ray slice is produced by `slices.one`, and the conformal slice by `View.slice`, so their moment versions come from the embedding data.

The curvature view reads the saved exact Ricci and Kretschmann scalars of the truncated metric through `curvature.nr_published`.
It holds $t=-\ell$ and plots them over $0.05\le r/\ell\le4$; it does not recompute curvature.
$K\ell^4$ has a logarithmic vertical scale from $10^{-4}$ to 1, and $R\ell^2$ a signed logarithmic scale from $-0.02$ to 1, linear within 0.01 of zero.
The shared curvature construction verifies the plotted samples at sixty digits.

No drawing kind is NOT_DRAWN for this page.
Only existing site classes are painted: `nr-world`; `cd-cover`, `cd-r`, `cd-t`, `cd-null`, `cd-centre`, `cd-chartedge` and `cd-scri`; the embedding's established cover, ring and meridian classes; and the curvature's standard K and R classes.
No glow is added.


## Physical checks

The printer's dedicated check compares the coefficient of $\epsilon$ in every Einstein-tensor slot with the O(4) sigma-model stress of $\chi=\pi+2\arctan(t/r)$, at three rational points and forty digits.
It separately checks the vacuum limit and the identically vanishing radial connection.
`test_global_texture.py` obtains the sphere area from the stored metric, derives $\Delta\Omega=4\pi(1-g_{\theta\theta}/r^2)$, and checks both limits.
Its growth proof uses $d[\arctan z-z/(1+z^2)]/dz=2z^2/(1+z^2)^2>0$ for $z>0$, with the bracket zero at $z=0$.
The test also checks the two source expressions at first order, the exact radial geodesics, and the spatial-volume factor underlying the test-dust density.

## Prose and browser checks

The collection's `view_text.length` counts are 193 characters for the radial-null view, 220 for the conformal view, 167 for the embedding, and 248 for the curvature view.
Each drawing kind has its own caption, and each aspect ratio lies between 0.5 and 2.
The headless local preview was checked at a 1200-pixel viewport and an emulated 390-pixel mobile viewport, without screenshots.
Both had no MathJax error elements and no document-wide horizontal overflow.
The conformal and embedding figures use their existing horizontal-scroll frames at the mobile width.
The three related links and their reciprocal links were followed successfully in the browser.
The painted curve classes resolve to existing styles, with no filter applied.
