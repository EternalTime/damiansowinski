# The torus universe in 2+1 dimensions

## The flat quotient

Steven Carlip's equation (1.7) in arXiv:hep-th/9109006v2 describes two commuting Poincaré transformations with rapidities $\lambda$ and $\mu$ and translations $a$ and $b$.
Let $x$ and $y$ have period one, let $t>0$, and choose the orientation $\Delta=a\mu-b\lambda>0$.
In units $c=1$, the developing map into Minkowski space is

$$T=t\cosh(\lambda x+\mu y),\qquad X=t\sinh(\lambda x+\mu y),\qquad Y=ax+by.$$

The induced spacetime metric is

$$ds^2=-dt^2+(a\,dx+b\,dy)^2+t^2(\lambda\,dx+\mu\,dy)^2.$$

The two spatial lattice vectors in the Euclidean development of a moment are $(a,t\lambda)$ and $(b,t\mu)$.
Their determinant is $t\Delta$, so the quotient is a compact flat torus with area $A=t\Delta$.
The complex modulus is

$$m(t)=\frac{b+it\mu}{a+it\lambda},\qquad \operatorname{Im}m=\frac{t\Delta}{a^2+t^2\lambda^2}>0.$$

These are Carlip's equations (1.9) and (1.11), with mean curvature $\tau=1/t$.
The period lattice is part of the global geometry.
A linear change of local spatial coordinates can diagonalize the metric while changing that lattice.

## Charts and sectors

The expanding sector has proper time $t>0$, the constant mean curvature coordinate $\tau=1/t$, and local Minkowski coordinates $(T,X,Y)$ with both deck transformations retained.
Time reversal supplies the contracting sector.
The separate static sector has two independent spatial translations and constant torus moduli.
Jorma Louko and Donald M. Marolf classify the two sectors with spacelike torus slices in section IV of arXiv:gr-qc/9308018v2.
Their equations (3.3) and (4.2) reproduce the evolving metric above, with $a_1=a$, $a_2=b$, $\lambda_1=\lambda$, and $\lambda_2=\mu$.
Their equations (3.6) and (4.1) give the static sector,

$$ds^2=-dt^2+\left(r_1\,dx+r_2\cos\phi\,dy\right)^2+r_2^2\sin^2\phi\,dy^2,$$

with $r_1,r_2>0$ and $0<\phi<\pi$.
Its two lattice vectors are $(r_1,0)$ and $(r_2\cos\phi,r_2\sin\phi)$.

Louko and Marolf's equation (4.3) supplies a further chart of the evolving sector:

$$ds^2=-du\left(\lambda\,dx+\mu\,dy\right)+u\left(\lambda\,dx+\mu\,dy\right)^2+\left(a\,dx+b\,dy\right)^2.$$

Here $u>0$ has dimensions of length squared.
For any length $\ell>0$, set $q=\lambda x+\mu y$ and

$$T=\frac{\ell e^q+u e^{-q}/\ell}{2},\qquad X=\frac{\ell e^q-u e^{-q}/\ell}{2},\qquad Y=ax+by.$$

Then $T^2-X^2=u$, and $-dT^2+dX^2=-du\,dq+u\,dq^2$.
The induced tori remain spacelike for $u>0$.
Continuation to $u<0$ makes the tori timelike, as Louko and Marolf show in section IV; the compact spatial cosmology occupies $u>0$.

## Intrinsic and spacetime curvature

The developing map is a local isometry to Minkowski space.
The spatial lattice has the Euclidean metric at each fixed positive $t$.

## Quotient drawings

The flat spatial torus has a parallelogram development with opposite sides identified.
A smooth doughnut surface in Euclidean 3-space has nonzero Gaussian curvature and cannot carry this flat metric globally.
Any spatial development must retain its paired edges.

In fact, no compact flat torus admits a smooth isometric immersion in Euclidean 3-space.
For a smooth compact surface immersed in Euclidean 3-space, choose a point of maximal Euclidean distance $R$ from an origin outside the surface.
At that point the Hessian of the squared distance along a unit tangent vector $v$ is $2+2R\,\mathrm{II}(v,v)\leq0$, with the normal directed radially outward.
Both principal curvatures therefore have the same nonzero sign, so the Gaussian curvature is positive there.
This contradicts the vanishing Gaussian curvature of a flat torus.
This obstruction is specific to Euclidean 3-space; a rectangular flat torus embeds smoothly and isometrically in Euclidean 4-space as a product of circles.
A planar fundamental parallelogram with paired edges is an intrinsic development.

A projection to one spatial coordinate suppresses one torus cycle and its null winding classes.
The generic period lattice supplies no preferred radial reduction retaining both cycles.
The two spatial identifications belong in a spacetime drawing with both spatial coordinates.

## Printed coordinates and figure construction

The printed proper time and static charts keep $c$ explicit and use tensor coordinates $x^0=ct$.
The mean curvature coordinate is $	au=1/(ct)$ and has dimensions of inverse length.
The null boundary coordinate $u$ has dimensions of length squared and is not scaled by $c$.
The five numerical chart checks compare pullbacks of the source's developing maps at three distinct points each.
No symbolic determinant or symbolic simplification enters those checks.

The full evolving family retains the four constants $a,b,\lambda,\mu$, with $a\mu-b\lambda>0$.
Neither rapidity nor either transverse translation is set to zero in the printed family.
The static family retains both side lengths and their angle.
An integer change of spatial basis of determinant one acts on the lattice and on the modulus by a modular transformation, while preserving the area and the quotient geometry.
The dedicated tests reconstruct the area and the complex modulus from the printed spatial Gram matrix, and hold its trajectory to a hyperbolic geodesic.
They also test the two commuting deck generators and closed timelike integer winding curves in the extension $u<0$.

Every chart has one three dimensional spacetime figure, with both pairs of opposite spatial cuts retained.
The evolving figures use $a=\ell$, $b=\ell/3$, $\lambda=1$, and $\mu=\sqrt{2}$, from $ct/\ell=1/2$ to $2$.
The static figure uses $r_1=\ell$, $r_2=11\ell/10$, and $\phi=\pi/3$.
The two edge colours name $x=0\sim1$ and $y=0\sim1$.
In each spatial development, the lattice vectors are drawn with their full length and angle at that moment.
The Minkowski figure uses the actual developing map on a centred cell, so its identified edges are related by the boosts and transverse translations.
The other four figures place intrinsic spatial developments at a height proportional to their displayed clock coordinate.
Every null cone is transformed by the Jacobian of that drawing map and checked against the full published metric and inverse metric.
The mean curvature figure takes decreasing $\tau$ as future; the null boundary figure uses the future normal to $u$.
The generated line coordinates are checked for both identification pairs in every chart.
The display windows are enlarged to keep every figure between 1:2 and 2:1.
These are spacetime figures with a continuous time axis, with no sequence of played spatial moments.

`embedding.NOT_DRAWN` records the obstruction proved above to a smooth global isometric immersion in Euclidean three-space.
The planar paired-edge developments in the spacetime figures are intrinsic cut descriptions of the flat spatial quotient.
`conformal.NOT_DRAWN` records the lack of a two-coordinate radial reduction retaining both independent winding classes.
This concerns a global two dimensional conformal diagram of the generic torus cosmology; local conformal flatness still holds.
The two curvature views are straight zero curves read from the stored Ricci and Kretschmann scalars through `curvature.py`.
No curvature is recomputed for either figure.

## Source access

The original full texts of Moncrief, Hosoya and Nakao, and Witten remain unread source debts.
Their verified primary records and the readable primary treatments by Carlip and by Louko and Marolf support the historical attributions and the classical chart family, as recorded in the task's literature ledger.
The detailed classical formulae above are taken from the fully read accessible treatments.

## Sources

Emil Martinec, Soluble systems in quantum gravity, Physical Review D 30, 1198-1204 (1984), section 5.
Vincent Moncrief, Reduction of the Einstein equations in 2+1 dimensions to a Hamiltonian system over Teichmüller space, Journal of Mathematical Physics 30, 2907-2914 (1989).
Steven Carlip, (2+1)-Dimensional Chern-Simons Gravity as a Dirac Square Root, Physical Review D 45, 3584-3590 (1992), arXiv:hep-th/9109006v2, section 1.
Jorma Louko and Donald M. Marolf, The solution space of 2 + 1 gravity on $\mathbb{R}\times T^2$ in Witten's connection formulation, Classical and Quantum Gravity 11, 311-330 (1994), arXiv:gr-qc/9308018v2, sections III and IV.
Steven Carlip, Quantum gravity in 2+1 dimensions: The case of a closed universe, Living Reviews in Relativity 8, 1 (2005), sections 2.2-2.8.
