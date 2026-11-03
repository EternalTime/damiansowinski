# Black hole with skyrmion hair

## Sources and chart

The static spherical exterior is written in the areal radius by Shiiki and Sawado, arXiv:gr-qc/0501025, equations (21)-(35).
The metric uses $c=1$ with $G$ explicit.

## Boundary value problem

Write $x=aF_\pi r$, $\mu=aF_\pi Gm$, $\alpha=4\pi GF_\pi^2$, and $C=1-2\mu/x$.
A prime below means $d/dx$.
Shiiki and Sawado's equations (30) and (31) become

$$\frac{N'}N=\frac{\alpha}{4x}(x^2+8\sin^2 f)f'^2,\qquad \mu'=\frac{\alpha}{8}\left[C(x^2+8\sin^2 f)f'^2+2\sin^2 f+\frac{4\sin^4 f}{x^2}\right].$$

The static field equation (32) is equivalently $(NCu f')'=N\sin 2f(1+4Cf'^2+4\sin^2 f/x^2)$ with $u=x^2+8\sin^2 f$.
The horizon has $\mu(x_h)=x_h/2$, and equations (33) and (34) give its first derivatives in `horizon_data`.
At infinity $N=1$ and $f=0$; the massless field has the decaying tail $f\sim x^{-2}$, as Dvali and Gußmann state in section II.C.

The example uses $\alpha=0.01$ and $x_h=0.1$ and shoots on the branch with larger $f_h$.
The integration begins at $x_h+10^{-6}$ with the horizon's first-order expansion and imposes $xf'+2f=0$ at $x=60$; the lapse is normalized there.
A convergence comparison at $x=30,60,120$ and offsets $10^{-5},10^{-6},10^{-7}$ gives $f_h=3.03711134,3.03711116,3.03711115$ for offset $10^{-6}$, respectively, with finite-radius masses $0.07748521,0.07748926,0.07748976$.
These are computed diagnostics, not numbers taken from a paper, and the finite outer boundary leaves a numerical approximation.
`jet` differentiates the same equations to provide the first and second derivatives of the numerical metric functions for the figures.

## Drawings

The areal chart is the exterior $r>r_h$ of the spherical solution.
The conformal diagram is not drawn: the sourced numerical chart ends at the horizon and does not establish the causal structure inside it.

## Checks
