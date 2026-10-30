# The charged rotating bubble

## The matching surface

Carlos A. López's 1984 model joins a flat interior to a Kerr-Newman exterior across a timelike oblate shell.
The readable treatment is Burinskii, Elizalde, Hildebrandt, and Magli, Physical Review D 65, 064039 (2002), arXiv:gr-qc/0109085v2, equations (2), (3), (16), (43), and (108).
Their mass length is $m = GM/c^2$, their signed charge length satisfies $q^2 = GQ^2/(4\pi\epsilon_0c^4)$, and $a = J/(Mc)$.
We use $m > 0$, $q^2 > m^2$, and $a \ge 0$.
The charge sign enters the electromagnetic potential and leaves the geometry invariant.

With $\Sigma = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 + a^2 - 2f$, the exterior has $f = mr - q^2/2$ and the interior has $f = 0$.
Continuity fixes $r_0 = q^2/(2m)$.
In the auxiliary Cartesian coordinates, the shell satisfies

$$\frac{x^2+y^2}{r_0^2+a^2}+\frac{z^2}{r_0^2}=1.$$

Its equatorial semiaxis is $\sqrt{r_0^2+a^2}$ and its polar semiaxis is $r_0$.
The full polar diameter $2r_0 = Q^2/(4\pi\epsilon_0Mc^2)$ is the classical electron radius when $M$ and $Q$ take electron values.
The oblate coordinate degeneracy at $r=0$ belongs to flat space inside the shell; Cartesian coordinates cover that interior without retaining a Kerr ring.

## Charts

Boyer-Lindquist coordinates follow equation (108), with $r>r_0$ outside and $0<r<r_0$ inside.
Kerr angular coordinates follow equation (94), with the same two radial regions and the transformation in Appendix C.
The flat interior also has Cartesian coordinates.
The shell is a separate distributional surface, excluded from the smooth coordinate domains.

## Surface stresses

The shell stress follows by taking the thin limit of equations (10)-(15), or by the jump in extrinsic curvature of equation (108).
Put $B_0=r_0^2+a^2$ and $\Sigma_0=r_0^2+a^2\cos^2\theta$.
The surface metric is $h_{ab}dx^adx^b=-(dx^0)^2+\Sigma_0d\theta^2+B_0\sin^2\theta d\phi^2$, with $x^0=ct$.
The outward normal on both sides points towards increasing $r$.
With $K_{ab}=\tfrac12n^r\partial_rg_{ab}$, exterior minus interior, the jump is

$$[K_{ab}]=\frac{m\sqrt{B_0}}{\Sigma_0^{3/2}}v_av_b,\qquad v_adx^a=dx^0-a\sin^2\theta d\phi,\qquad [K]=-\frac{m}{\sqrt{B_0\Sigma_0}}.$$

Israel's condition in the $(-,+,+,+)$ convention is $S_{ab}=-([K_{ab}]-h_{ab}[K])/(8\pi)$.
Consequently,

$$u^a=\frac{(B_0,0,a)}{\sqrt{B_0\Sigma_0}},\qquad S_{ab}=-\mathcal T(h_{ab}+u_au_b),\qquad\mathcal T=\frac{m}{8\pi\sqrt{B_0\Sigma_0}}.$$

Here $S_{ab}$ is the geometric surface stress; multiplying it by $c^4/G$ restores physical energy per area.
The comoving surface energy density is zero and both tangential pressures equal $-\mathcal T$.
These are derived from the readable 2002 metric, not transcribed from the unread 1984 or 1985 texts.

The independent thin limit uses $f=(mr-q^2/2)\Theta(r-r_0)$, so $f'=m\Theta(r-r_0)$ and $f''=m\delta(r-r_0)$.
Equations (10)-(15) then give no delta term in the comoving density and $-m\delta(r-r_0)/(8\pi\Sigma_0)$ in each tangential pressure.
Since $d\ell=\sqrt{\Sigma_0/B_0}\,dr$, the coefficient of $\delta(\ell)$ is precisely $-\mathcal T$.
The 2002 authors use the opposite Ricci convention and write Einstein's equation with $-R_{ab}$; the physical stress and its signs agree with the standard Ricci contraction used here.
Every smooth bulk tensor excludes the shell.
A product such as the square of its delta curvature is undefined as an ordinary distribution, so no shell value of the Kretschmann scalar is assigned.

The angular velocity is $d\phi/dt=ca/B_0$, and the speed measured in the interior inertial frame is $ca\sin\theta/\sqrt{B_0}$.
It is subluminal at every latitude.
The dedicated numerical check of the extrinsic curvature and the matching remains to be completed.

## Electromagnetic boundary

The exterior potential of Burinskii 2016, arXiv:1602.04215v1, equations (4)-(5), has intrinsic form $A^+=-qr\Sigma^{-1}(dx^0-a\sin^2\theta d\phi)$ in the time and spin orientation of the 2002 metric.
The 2016 time is reversed to this orientation; the Kerr angular and Boyer-Lindquist azimuths agree on $dr=0$.
For a field-free interior, direct differentiation on the shell gives

$$\iota^*[F]=\frac{2qr_0a\sin\theta\cos\theta}{\Sigma_0^2}\,d\theta\wedge(B_0d\phi-a\,dx^0).$$

This is nonzero at generic latitude for $qa\ne0$.
In the co-rotating azimuth $\varphi=\phi-(a/B_0)x^0$, it is magnetic flux through the material surface.
The total flux over the shell vanishes, and its local value is nonzero.
An ordinary electric surface current leaves the homogeneous Maxwell equation $dF=0$ unchanged, which requires the intrinsic jump to vanish for a bounded field with a discontinuity.
Thus a field-free interior and the rotating exterior require additional electromagnetic boundary structure.

A formal electric dipole layer can be introduced through $A=\Theta(r-r_0)A^+$, whose field contains $\delta(r-r_0)dr\wedge\iota^*A^+$.
The electromagnetic stress of that singular field contains delta products and needs a prescribed regularization.
For this explicit potential step, the quadratic invariant has coefficient $-2q^2r_0^2/\Sigma_0^2$ multiplying $\delta(r-r_0)^2$.
The gravitational surface tensor fixes the total stress required by the geometry and does not specify this electromagnetic prescription.
A magnetic surface source is another idealization; a nonzero interior Maxwell field would require compensating matter stress to retain flatness.
The separate 2016 model uses Higgs/chiral fields, a finite transition zone, and penetrating currents, equations (12), (17), and (30).

## Causal domain

The restriction $q^2>m^2$ makes $\Delta=(r-m)^2+a^2+q^2-m^2$ positive.
For $r\ge r_0$, $f\ge0$ and $g^{tt}<0$ in Boyer-Lindquist coordinates.
The common time therefore extends from the exterior to the flat interior with timelike gradient.
The physical model retains one exterior and a flat ball, with no negative radial sheet.

## Drawings

Parameters, null slices, embedding, and curvature views remain to be selected from this causal domain.
