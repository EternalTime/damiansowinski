# Universes of stiff fluid

## Scope and sources

The page publishes an independently derived plane-symmetric stiff-fluid family.
The verified founding records are Tabensky and Taub (1973), DOI 10.1007/BF01661153, and Wainwright, Ince and Marshman (1979), DOI 10.1007/BF00759860.
Their abstracts were read; their full texts were inaccessible through the legitimate routes recorded in the external task literature ledger.
Rendall (1994 preprint, 1995 journal), equations (2)-(5), supplies the areal reduction.
The coordinate transformations, stress, curvature, and drawings below are checked independently.

## Normalization

Use signature (-,+,+,+), G=c=1, and canonical scalar Phi with Einstein coupling eight pi explicit.
The dimensionless review scalar sigma is sqrt(8 pi) times Phi.
The scalar's stress equals a stiff fluid's stress wherever its gradient is timelike.

## Charts

With a fixed length unit $\ell$, all coordinates dimensionless, the areal line element is

$$ds^2=\ell^2[F(-dt^2+dx^2)+t(dy^2+dz^2)],\qquad F=t^{a^2-1/2}\exp(b^2t^2/2+2abx).$$

The review scalar is $\sigma=a\ln t+bx$ and the canonical field is $\Phi=\sigma/\sqrt{8\pi}$.
Direct differentiation gives $R_{tt}=a^2/t^2$, $R_{tx}=ab/t$, $R_{xx}=b^2$, with every other Ricci component zero.
The scalar wave equation reduces to $\Phi_{tt}+\Phi_t/t=\Phi_{xx}$.
The norm is $\nabla\sigma\cdot\nabla\sigma=-(a^2/t^2-b^2)/(\ell^2F)$, so the fluid has

$$\rho=p=\frac{a^2/t^2-b^2}{16\pi\ell^2F},\qquad U^\mu=-\frac{\nabla^\mu\sigma}{\sqrt{-\nabla\sigma\cdot\nabla\sigma}}.$$

Then $T_{\mu\nu}=2\rho U_\mu U_\nu+\rho g_{\mu\nu}=\Phi_{,\mu}\Phi_{,\nu}-g_{\mu\nu}(\nabla\Phi)^2/2$.
The fluid occupies $0<t<a/b$ for $b>0$ and all $t>0$ for $b=0$.
At $t=a/b$ the stress remains finite at finite $x$ even though the unit fluid velocity has no timelike limit.
The scalar spacetime continues smoothly across that boundary.
At $t=0$ the Ricci scalar diverges for all $a>0$.

Comoving labels are $\chi=x+bt^2/(2a)$, whose gradient annihilates the fluid velocity.
Substituting $dx=d\chi-(bt/a)dt$ changes $F$ to $t^{a^2-1/2}\exp(-b^2t^2/2+2ab\chi)$ and introduces the cross term $-2Fbt\,dt\,d\chi/a$.
These labels follow the fluid while retaining areal time; they are independently derived and are not claimed to reproduce an unread founding chart.
The characteristic coordinates are $u=t-x$, $v=t+x$, giving $ds^2=\ell^2[-Fdu\,dv+(u+v)(dy^2+dz^2)/2]$.
The null chart is shown against $(v-u)/2$ and $(u+v)/2$ in its ray diagram.

The conformal factor is held with explicit rational slopes while the tensors are printed.
The chart checks compare every Ricci slot and each slope at points, without a symbolic determinant or simplification inside a check.
The stored Kretschmann scalar in the areal and comoving charts is

$$K=\frac{3+8a^4-4a^2+8b^4t^4+4b^2t^2-16a^2b^2t^2}{4F^2\ell^4t^4}.$$

At $b=0$, proper time is $\tau=\ell t^{a^2/2+3/4}/(a^2/2+3/4)$.
The three powers of $\tau$ are $(a^2-1/2)/(a^2+3/2)$ and twice $1/(a^2+3/2)$.
They sum to one and all become $1/3$ at $a^2=3/2$.
The five independent tests read the published metric and Einstein tensor and test the canonical scalar, the fluid stress, the dictionary, the boundary, the coordinate changes and the isotropic powers.

## Drawings

The drawings use a=b=1 and the fluid region 0<t<1.
The scalar extension beyond its upper boundary is regular at finite $x$.
The ray views span $|x|\leq0.8$ (or $|\chi|\leq0.8$) and $-0.3\leq t\leq1.3$, including generous margins beyond the fluid strip.
A marked fluid world line has $\chi=0$.
The conformal map is $p=\arctan(4(t-x))$, $q=\arctan(4(t+x))$, with $X=q-p$ and $T=q+p$.
It maps the strip onto the lens above the singular line $t=0$ and below the regular scalar boundary $t=1$.
Every chart's map is held against its own published metric and future timelike direction.
The embedding is one symmetry plane at $t=1/2$, $x=0$, with induced metric $\ell^2t(dy^2+dz^2)$.
Its disc has proper radius $\ell$ and is intrinsically flat; the generator checks the scale $\ell\sqrt{t}$ and every profile chord.
All these views are single frames, so no movie or stack entry is needed.

## Verification

Every chart is checked at points against the scalar and perfect-fluid stresses, the field dictionary, and the fluid boundary.

## Curvature maps

Both curvature maps use the stored scalars of the areal chart at $a=b=\ell=1$, on $-0.5\leq x\leq0.5$, $0.05\leq t\leq1$.
The $80\times76$ cell values are evaluated from those expressions, with a sample checked again to sixty digits.
The $K$ colour range is $10^{-2}$ to $10^8$ and the $R$ range is $-10^4$ to zero, signed logarithmic with linear width $10^{-3}$.
Every sampled cell fits within its colour scale.
At fixed $t$, $K$ varies as $e^{-4x}$ and $R$ as $e^{-2x}$, checked from the stored scalars.
As $t\to0$, their magnitudes diverge as $t^{-5}$ and $t^{-5/2}$, respectively.
At $t=1$, $R=0$ and $K=3e^{-1-4x}/4$, finite and positive at finite $x$.
The embedding plane meets each quotient diagram at one point, $(t,x)=(1/2,0)$, transformed into the other charts by their own dictionaries.

The characteristic domain explicitly includes $v>-u$ so that $t=(u+v)/2>0$.
That dependent lower bound also keeps automatic symmetry sampling inside the Lorentzian patch, where the transverse Killing vectors have positive norm.
The remaining upper fluid bound is the declared inequality $u+v<2a/b$ when $b>0$.
