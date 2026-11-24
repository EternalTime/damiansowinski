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

## Embedding diagram

### What the metric has

A moment of $t$ is the same surface in the areal and comoving charts, and in the comoving chart its metric is $\ell^2[F\,d\chi^2+t(dy^2+dz^2)]$ with $F=t^{a^2-1/2}e^{2ab\chi-b^2t^2/2}$.
Every moment is flat.
The proper distance from the fluid world line $\chi=0$,

$$X=A(t)\,\frac{e^{ab\chi}-1}{ab},\qquad A(t)=\ell\,t^{a^2/2-1/4}e^{-b^2t^2/4},$$

turns it into $dX^2+\ell^2t(dy^2+dz^2)$, Euclid's space with its transverse axes scaled by $\ell\sqrt t$.
The factor $e^{2ab\chi}$ is the whole of the dependence on $x$, and $X$ absorbs it, so the geometry of a moment varies neither along $x$ nor within the planes.
The variation along $x$ lies in the density, $\rho\propto e^{-2ab\chi}$, and in the curvature, $K\propto e^{-4ab\chi}$, which the curvature maps draw.
Each moment ends where $\chi\to-\infty$, at the finite proper distance $A(t)/(ab)$ from $\chi=0$, where the density and the curvature diverge; that plane lies at infinite coordinate distance and no causal curve inside the fluid region reaches it.
The History cites neither fact, so no caption states them.

The fluid moving through the moments shows the dependence on $t$.
The fluid elements keep their labels $\chi$, $y$, and $z$: the raised gradient of $\sigma=a\ln t+b\chi-b^2t^2/(2a)$ has no $\chi$, $y$, or $z$ component.
In a moment of $t$ the proper distance between two fluid elements across the planes is $A(t)$ times a number fixed by their labels, and within the planes $\ell\sqrt t$ times one.
So a ring of the fluid round at one moment is at every other moment an ellipse whose semi-axes scale as $A(t)$ across the planes and as $\sqrt t$ within them.
A shift of $\chi$ multiplies the whole quotient metric by a constant, so the shapes the ring takes are the same about every fluid world line.
For $a^2>1/2$ the ring widens across the planes until $t=\sqrt{a^2-1/2}/b$ and narrows afterwards, toward $t=a/b$, where the fluid crosses the planes at light speed.
At $b=0$ the same ring widens across the planes as $t^{a^2/2-1/4}$ for all $t$, as in Kasner's universe with a scalar field and equal transverse exponents.

### The decision

The embedding diagram is a ring of the fluid in the plane $z=0$, in two views, the form Kasner's universe with a scalar field takes: the ring's world tube, a stack of ellipses, and a movie of the flat moments with the ring on them.
The ring is of the fluid, the world lines of constant $\chi$, $y$, and $z$ of the comoving chart, so its stretch is the solution's own pair of scale factors.
The fluid is not free: its four acceleration, computed from the comoving metric, is $b$ times a factor that at $a=b=1$ vanishes nowhere in $0<t<1$, so the ring is named a ring of the fluid.
A ring of free particles would need a moment at which it starts at rest, and the metric fixes none; particles at rest at $t=1/2$ were crossing the planes at $0.75$ of light speed at $t=0.01$.
The transverse plane alone is the flat plane $\ell^2t(dy^2+dz^2)$ at every $x$, the picture of flat space, and a movie of it shows a disc growing as $\sqrt t$ and nothing of this spacetime's own.
No quantity is drawn as a height, since the density and the curvature already have their maps.

The parameters are $a=b=\ell=1$, the values of every other drawing of this spacetime.
The movie runs from $t=0.1$ to $0.95$ inside the fluid region $0<t<1$, a frame every $0.025$, 35 frames, with the stack's five moments $t=0.1$, $0.25$, $0.5$, $0.75$, and $0.95$.
At $t=0.95$ the fluid moves across the planes at $0.95$ of light speed in the areal chart, and $t=1/\sqrt2$, where the ring is widest across the planes, lies between the moments $0.5$ and $0.75$.
The ring is the circle of radius $\ell/3$ about $\chi=0$ at $t=1/2$; at $t=0.1$ its semi-axes are $0.24\,\ell$ across the planes and $0.15\,\ell$ within them, and at $t=0.95$ they are $0.33\,\ell$ and $0.46\,\ell$.
Each moment is drawn as a disc of radius $\ell/2$ about the ring's centre, smaller than $A(t)/(ab)\geq0.56\,\ell$ at every frame, so every point drawn belongs to the moment.
Twelve fluid elements are marked, and their world lines run up the tube, which stands $1.5\,\ell$ high for each unit of $t$ on the axis $\chi=0$.

### What the generator checks

`StretchedPlane` in `embedding.py` reads the plane of $\chi$ and $y$ from the published comoving metric, with the chart's own $F$ written out, and takes $X(\chi)$ as the integral of $\sqrt{g_{\chi\chi}}$ in sympy, checked by differentiating it.
Each frame is held to the published metric along, across, and around, as every slice is, and to the closed form $A(t)/(ab)$ for the distance to $\chi\to-\infty$, which must exceed the disc's radius.
The ring's labels are found from its circle at $t=1/2$, and at every frame the ring drawn from the published metric is checked against the quadrature of $\sqrt{g_{\chi\chi}}$ from $\chi=0$ and against the ellipse of semi-axes $(\ell/3)(2t)^{1/4}e^{(1-4t^2)/16}$ and $(\ell/3)\sqrt{2t}$.
`stiff_fluid_flow` checks the published Einstein tensor against the stress of $\sigma$, the raised gradient of $\sigma$ for its zero $\chi$, $y$, and $z$ components, and its norm for being timelike at the moments drawn.
Every row of the tube is checked against the ring drawn from the published metric at its own time.

### The moments on the other drawings

Each of the five moments is the line of constant $t$ across the disc it is drawn as.
The profile ends at $\chi_{hi}$ with $e^{\chi_{hi}}-1=\tfrac12/A(t)$, and the disc's far side lies at $e^{\chi_{lo}}=2-e^{\chi_{hi}}$.
That reach is $\chi$ in the comoving chart, $x=\chi-t^2/2$ in the areal chart, and $(u,v)=(t-x,t+x)$ in the characteristic chart, and on the conformal diagram the curve $p,q=\arctan(4(t\mp x))$ over the same $x$.
`slices.py` checks both maps by pulling the areal metric back onto the comoving and characteristic ones, and every profile against $X=t^{1/4}e^{-t^2/4}(e^\chi-1)$.
A moment of the tube is the moment of the movie at its time, and both views mark all five.

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

The characteristic domain explicitly includes $v>-u$ so that $t=(u+v)/2>0$.
That dependent lower bound also keeps automatic symmetry sampling inside the Lorentzian patch, where the transverse Killing vectors have positive norm.
The remaining upper fluid bound is the declared inequality $u+v<2a/b$ when $b>0$.
