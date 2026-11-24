# Universes with a scalar field and an exponential potential

## Sources and normalization

We take $c=1$, $\kappa^2=8\pi G$, and a canonical minimally coupled scalar with $V=V_0e^{-\lambda\kappa\phi}$, $V_0>0$, and $\lambda>0$.
Copeland, Liddle and Wands (1998), equations (1)-(4), Table I and (A7)-(A8), and Heard and Wands (2002), equations (1)-(3), (10)-(11) and (22)-(24), fix the matter normalization and the two exact branches.
The equations are $3H^2=\kappa^2(\rho_\gamma+\dot\phi^2/2+V)$, $\dot H=-\kappa^2(\gamma\rho_\gamma+\dot\phi^2)/2$, $\ddot\phi+3H\dot\phi-\lambda\kappa V=0$, and $\dot\rho_\gamma=-3\gamma H\rho_\gamma$.
The uncoupled fluid has pressure $p_\gamma=(\gamma-1)\rho_\gamma$.

## Exact matter branches

With $a=a_0(t/t_0)^p$ and $\phi=2\ln(t/t_0)/(\lambda\kappa)$, the scalar equation fixes $V_0t_0^2=2(3p-1)/(\kappa^2\lambda^2)$.
The scalar density is $\rho_\phi=6p/(\kappa^2\lambda^2t^2)$ and the fluid density is $\rho_\gamma=3p(p-2/\lambda^2)/(\kappa^2t^2)$.
Positive potential requires $p>1/3$.
The zero of the scalar is fixed at $t_0$, tying the amplitude to that time.

The scalar only branch has $p=2/\lambda^2$, $0<\lambda^2<6$, and $V_0t_0^2=2(6-\lambda^2)/(\kappa^2\lambda^4)$.
Its equation of state is $w_\phi=\lambda^2/3-1$ and it accelerates for $\lambda^2<2$.
Homogeneous perturbations in the expanding flat scalar only model have the eigenvalue $(\lambda^2-6)/2$ per logarithmic scale factor.
With a fluid present, the extra eigenvalue is $\lambda^2-3\gamma$, so attraction in that model requires $\lambda^2<3\gamma$.

The fluid scaling branch has $0<\gamma<2$, $p=2/(3\gamma)$, and $\lambda^2>3\gamma$.
Its amplitude is $V_0t_0^2=2(2-\gamma)/(\gamma\kappa^2\lambda^2)$ and its fluid density is $4(\lambda^2-3\gamma)/(3\gamma^2\kappa^2\lambda^2t^2)$.
The scalar fraction is $\Omega_\phi=3\gamma/\lambda^2$ and $w_\phi=\gamma-1$.
Both eigenvalues in Copeland, Liddle and Wands (A8) have negative real part in this open range.
At $\lambda^2=3\gamma$ the fluid vanishes and linear stability is degenerate.
Acceleration requires $\gamma<2/3$.

Wetterich (1988), (53)-(59) and (64), uses $\phi=2\sqrt{3+2\omega}\,\sigma$, $M^2=1/(2\kappa^2)$, and $\lambda^2=A^2/[2(3+2\omega)]=2/Y$.
His $n=3\gamma$ and zero direct matter coupling identify his fraction $nY/2$ with $3\gamma/\lambda^2$.

## Coordinate charts and causal domains

The cosmic Cartesian chart is the flat Friedmann geometry used by Copeland, Liddle and Wands and Heard and Wands.
The cosmic spherical chart follows by $x=r\sin\theta\cos\phi$, $y=r\sin\theta\sin\phi$, and $z=r\cos\theta$.
Abbott and Wise (1984), (2.1)-(2.3), use the cosmic and conformal Cartesian charts for $p>1$.
Integrating $d\eta=dt/a$ yields $\eta=t_0(t/t_0)^{1-p}/[a_0(1-p)]$ for $p\ne1$.
The same integral is positive for $1/3<p<1$ and has range $(0,\infty)$.
For $p>1$ its range is $(-\infty,0)$.
Coasting has $p=1$, $\eta=(t_0/a_0)\ln(t/t_0)$, and $\eta\in\mathbb R$.

The particle horizon has proper radius $t/(1-p)$ for $p<1$.
The event horizon has proper radius $t/(p-1)$ for $p>1$.
Both defining integrals diverge for $p=1$.
These domains assume the exact expansion over $0<t<\infty$.
The spherical origin is a regular spatial axis at every positive time.
The curvature scalars are $R=6p(2p-1)/t^2$ and $K=12p^2[(p-1)^2+p^2]/t^4$.
Thus $t=0$ is a curvature singularity throughout $p>1/3$, including $p=1/2$ where $R=0$.
The strictly decreasing $K$ has a timelike gradient, so a Killing field must be tangent to a spacelike constant time slice.

## Drawings and checks

We choose $a_0=t_0=1$ for the drawings.
The accelerating scalar example has $p=2$, $\lambda=1$.
The dust scaling example has $p=2/3$, $\gamma=1$, $\lambda=3$, and $\Omega_\phi=1/3$.
The coasting scalar example has $p=1$, $\lambda^2=2$.
The conformal maps compactify $u=\eta-r$ and $v=\eta+r$ with $U=\arctan u$ and $V=\arctan v$.
Each point in a radial conformal diagram represents a 2-sphere, with $r=0$ the axis.
The finite conformal boundary is spacelike, future infinity for $p>1$ and the initial singularity for $p<1$.
For coasting the singularity lies at past null infinity and the past timelike corner of the compactified domain.

## Embedding diagram

### What the metric has

A moment of cosmic time is $a^2(dx^2+dy^2+dz^2)$ in the cosmic Cartesian chart, with $a=a_0(t/t_0)^p$, so every moment of every member is flat.
Its plane $z=0$ is Euclid's plane with every length scaled by $a$, and one moment alone is the picture of empty flat space for every member.
The members differ in how $a$ grows with $t$.
The published Christoffel symbols have no $\Gamma^i{}_{tt}$, so a free particle at rest in the chart stays at rest, and a ring of them on $x^2+y^2=R^2$ is at every moment the circle of proper radius $aR$.
At $a_0=t_0=1$ that radius grows as $t^2$ at $p=2$, faster and faster, as $t^{2/3}$ at $p=2/3$, more and more slowly, and as $t$ at $p=1$, at a steady rate.
The members differ in their horizons too, which the spacetime diagrams and the conformal diagram already mark.
For $p>1$ the event horizon of the observer at $r=0$ has the comoving radius $\int_t^\infty dt'/a=t^{1-p}/(p-1)$ and the proper radius $t/(p-1)$.
For $p<1$ the particle horizon has the comoving radius $\int_0^t dt'/a=t^{1-p}/(1-p)$ and the proper radius $t/(1-p)$.
At $p=1$ both integrals diverge, and the observer has neither horizon.

### The decision

The embedding diagram has three views, one for each member drawn on every other drawing, $p=2$, $2/3$, and $1$ at $a_0=t_0=1$.
Each view is a movie of the plane $z=0$ of the cosmic Cartesian chart through a ring of free particles at rest in the chart, the form of the Big Rip's ring of galaxies.
Every movie runs over the same span, $t_0/2\le t\le3t_0/2$, a frame every $t_0/40$, 41 frames, with the three key moments $t=t_0/2$, $t_0$, and $3t_0/2$.
Over that span the ring grows by a factor of nine at $p=2$, two at $p=2/3$, and three at $p=1$, and its growth speeds up, slows down, or holds steady from frame to frame, so the three movies differ on sight.
The accelerating and the decelerating members also show the horizon of the observer at the ring's centre: the event horizon at $p=2$, of proper radius $t$, and the particle horizon at $p=2/3$, of proper radius $3t$.
Each of those rings has the comoving radius at which it meets its horizon at $t=t_0$, the middle of the span: $R=1$ at $p=2$ and $R=3$ at $p=2/3$.
At $p=2$ the ring lies inside the event horizon before $t_0$ and outside it after; at $p=2/3$ it lies outside the particle horizon before $t_0$ and inside it after.
The coasting ring has $R=1$, of proper radius $t_0$ at $t_0$ as the accelerating ring is.
Each moment is the disc of proper radius $2.5\,t_0$, $5\,t_0$, and $1.75\,t_0$ about the ring's centre, which holds the ring and the horizon at every frame.
The horizon is a circle about the disc's centre, so it is drawn as a ring of the disc of class `horizon`, as Oppenheimer and Snyder's collapse marks its horizon on every frame, and each frame holds the disc's profile, the ring of particles, twelve particles marked every $30°$, and that one circle.
No world tube stands beside the rings: three tubes and three movies would make six views, and each movie shows its member's rate of growth.
No quantity is drawn as a height, since every moment is homogeneous and the curvature scalars have their own figures.
The captions state what the frames show: the span, each ring's growth, each horizon's proper radius, and the moment the two meet.

### What the generator checks

`embedding_views` in `exponential_scalar_universes.py` reads each moment as a `FlatPlane` from the published cosmic Cartesian metric, with the chart's own definition of $a$ written out, and checks its scale against $t^p$.
Each frame is held to the published metric along, across, and around, as every slice is, and its ring to the circle of proper radius $Rt^p$.
`horizon_radius` integrates $dt/\sqrt{g_{xx}}$ of the published metric in sympy toward the future and toward the past, and checks that at $p=2$ the first converges, to $1/t$, and the second diverges, that at $p=2/3$ the second converges, to $3t^{1/3}$, and the first diverges, and that both diverge at $p=1$.
Every horizon circle is checked against the proper radius $t/|p-1|$ and to lie on the disc.
The ring's radius at the frames is checked to grow faster and faster at $p=2$, more and more slowly at $p=2/3$, and steadily at $p=1$, and to meet each horizon at $t=t_0$, inside or outside it on either side as above.
`at_rest` checks the published Christoffel symbols for the missing $\Gamma^i{}_{tt}$.

### The moments on the other drawings

Each key moment is the whole plane $z=0$, so on every spacetime diagram of its member it is the line of constant time across the drawing: $t$ in the cosmic charts, $\eta=t^{1-p}/(1-p)$ in the conformal chart, and $\eta=\ln t$ in the coasting chart.
On each conformal diagram it is the curve of that $\eta$ over every $r\geq0$, $U,V=\arctan(\eta\mp r)$.
The three members are three spacetimes of one line element, so each drawing marks only its own member's moments, which `hidden_views` in its exceptions file records.
