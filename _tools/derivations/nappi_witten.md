# The closed universe of Nappi and Witten

## Sources and normalisation

The cosmological metric is equation (16) of hep-th/9206078.
We restore the common length scale as $L^2$ and use $b^2=(1-\sin\alpha)/(1+\sin\alpha)>0$.
The dilaton is $\Phi=\Phi_0-\tfrac12\log D$, with a constant absorbed into $\Phi_0$.
This is the limit of large levels of the gauged Wess-Zumino-Witten model.

## Coordinate transformations

The regular gauge of equation (21) has $x=\tan^2\psi$ and $\lambda'=\lambda-\log(x)/(2b)$ on the cosmological overlap.
The invariant variable of Bars and Sfetsos, hep-th/9208001 (18)-(23), is $r=\sin^2\psi=x/(1+x)$.
The static continuation is $r=-\sinh^2\chi$, the continuation stated on the last page of the 1992 discussion and written explicitly for $b=1$ in hep-th/0204189 (3.10).
The finite level bosonic metric in hep-th/9208001 (19)-(22) has different curvature and is a quantum correction to this classical geometry.

## Spatial volume and the turnaround

Set $z=b\tan\psi$ and $y=\log z$.
For angular periods $2\pi$, the spatial volume is $V=2\pi^2L^3y/\sinh y$, with its continuous value $2\pi^2L^3$ at $y=0$.
For $b=1$, the maximum is at $\psi=\pi/4$, where $R=14\cos^2(2s)/L^2$.
The maximum for arbitrary $b$ is at $\psi_*=\arctan(1/b)$.
Since $y'=1/(\sin\psi\cos\psi)>0$, the sign of $V'$ is the sign of $\sinh y-y\cosh y$.
Its derivative with respect to $y$ is $-y\sinh y<0$ away from zero, so the volume increases up to $y=0$ and decreases after it.
Both endpoints have zero limiting volume.
At the maximum, $V''=-(2\pi^2L^3/3)(b+b^{-1})^2<0$.

For $b=1$, write $C=\cos(2s)$ at $\psi=\pi/4$.
The spacetime curvature scalars are $R=14C^2/L^2$ and $K=(76C^4+128C^2-128)/L^4$.
They are finite on the whole spatial slice, including the regular limits at its circular axes.
The Kretschmann scalar is negative near $s=\pi/4$; in Lorentzian signature it is a contraction with signs and can be negative.
The two circle expansion rates are $H_\rho=-(1+C)/L$ and $H_\lambda=(1-C)/L$.
Thus the total volume maximum coexists with local expansion $\Theta=-2C/L$ of either sign.

## Background fields

For the cosmological coordinates, a gauge potential is $B_{\rho\lambda}=L^2\sin^2s\cos^2\psi/D$.
It differs from the potential of equation (14) by a constant antisymmetric tensor.
Its field strength $H=dB$ satisfies $R_{\mu\nu}+2\nabla_\mu\nabla_\nu\Phi-\tfrac14H_{\mu\alpha\beta}H_\nu{}^{\alpha\beta}=0$ with the dilaton of equation (17).
The regular angle is denoted $w=\lambda'$.
In that gauge $B_{\rho w}=L^2\sin^2s/D$ and $B_{\rho x}=-L^2b\cos^2s/(2D)$.
The latter follows after subtracting a closed two form from the direct coordinate pullback, making the potential regular at $x=0$.

## The later plane wave

Equation (17) of hep-th/9310112 is $da_1^2+da_2^2+2du\,dv+(a_2da_1-a_1da_2)du+b_0du^2$.
The change $V=-v-b_0u/2$ makes it the rotating Cahen-Wallach form at $\omega=-1/2$ in the source's dimensionless units.
A rotation $X+iY=e^{-iu/2}(a_1+ia_2)$ makes the two transverse eigenvalues $\lambda_1=\lambda_2=-1/4$.
Its Ricci component is $R_{uu}=1/2$, in agreement with the source.
The standalone comparison `_tools/derivations/nappi_witten_wave.py` pulls the live Cahen-Wallach rotating metric back to the 1993 line element exactly.

## Drawing geometry

At $b=1$ and $\psi=\pi/4$, a surface of constant $\rho$ has metric $L^2(ds^2+\sin^2s\,d\lambda^2)$.
Its Euclidean surface of revolution has radius $L\sin s$ and height $L\cos s$, for $0\leq s\leq\pi/2$.
The same moment has $x=1$ in the regular gauge and $r=1/2$ in the invariant coordinate.
Each of the three cosmological null views marks the entire interval of $s$ reached by that hemisphere and names its embedding view `turnaround`.

The curvature curves use the stored four dimensional scalars on that same moment, so their horizontal coordinate is the polar angle $s$.
The constant $s$ plane of periodic time and $\chi$ in the static extension carries null curves; freely propagating light also moves across $s$.
The regular gauge's cosmological null view divides out $w$ and draws null geodesics of zero angular momentum.

NOT_DRAWN: a global conformal diagram through the static extensions would identify timelike separated points because their time is periodic, and the coset is non-Hausdorff at the joins.
The conformal construction is the closed cosmological causal region with both circular directions suppressed.
NOT_DRAWN: the maximum volume embedding moment has $r=1/2$ and lies on none of the static extension, where $r<0$.

The conformal null coordinates are $p=(\psi-s)/2$ and $q=(\psi+s)/2$, so the plotted axes are $X=q-p=s$ and $T=q+p=\psi$.

The null integrals in the transformed charts contain $\sqrt{x}$ and $\arcsin\sqrt{r}$.
The ray comparison omits $x\leq0.02$ and $r\notin(0.02,0.98)$, where these coordinate maps have unbounded derivatives at the endpoint.
Their tracing step is $0.001$ of the displayed box.
The static null integral is compared only for $\chi>0.02$, inside its sourced domain.
