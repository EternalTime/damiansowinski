# Extended inflation

## Sources and normalisation

Kolitch and Eardley, Annals of Physics 241 (1995), equations (2.2)-(2.8) and (5.20)-(5.21), supply the readable primary equations.
Their exact flat solution is $a=a_0(1+\chi t_s)^p$, $\varphi=\varphi_0(1+\chi t_s)^2$, $p=\omega+1/2$, and $\chi^2=32\pi\rho_v/[\varphi_0(6\omega+5)(2\omega+3)]$ in units $c=1$.
Set $t=t_s+1/\chi$, $t_0=1/\chi$, and absorb $a_0$ into the comoving radius.
With $\varphi_0=1/G_0$, the dimensionless scalar is $\Phi=(t/t_0)^2$, while $a=(t/t_0)^p$.
Restoring $c$ for an energy density gives $\rho_v=c^2(6\omega+5)(2\omega+3)/(32\pi G_0t_0^2)$.
The associated constant $\Lambda_v=8\pi G_0\rho_v/c^4=(6\omega+5)(2\omega+3)/(4c^2t_0^2)$ has dimensions of inverse length squared.
The Jordan action, in $x^0=ct$ coordinates, is $S=c^3/(16\pi G_0)\int\sqrt{-g}[\Phi R-\omega(\nabla\Phi)^2/\Phi-2\Lambda_v]\,d^4x$.
The constant false vacuum term corresponds to $T_{\mu\nu}=-\rho_vg_{\mu\nu}$.
The exact source normalization comes from this later primary derivation; the 1984 paper, 1989 Letter and its erratum, and Weinberg's 1989 full text remain unread access debts recorded in the task literature ledger.

## Charts

The cosmic spherical chart has $ds^2=-c^2dt^2+a^2(dr^2+r^2d\Omega^2)$, $t>0$, $t_0>0$, and $p>1$.
Both $t$ and $t_0$ are times, $r$ is a length, and $a$ is dimensionless.
Integrating Kolitch and Eardley's section III definition $d\eta=dt/a$ gives $\eta=-t_0(t/t_0)^{1-p}/(p-1)<0$.
Thus $a=[-(p-1)\eta/t_0]^{-p/(p-1)}$ and $ds^2=a^2(-c^2d\eta^2+dr^2+r^2d\Omega^2)$.
Both charts cover the same homogeneous branch.
Their components use $x^0=ct$ and $x^0=c\eta$, respectively; the generator's internal scale definitions include that change of units.
The exact powers are held as $a$ while their first derivatives, $\partial_{ct}a=pa/(ct)$ and $\partial_{c\eta}a=-pa/[(p-1)c\eta]$, determine the tensors.

## Field equations

The equations checked are $G_{\mu\nu}=-(\Lambda_v/\Phi)g_{\mu\nu}+\omega[\Phi_{,\mu}\Phi_{,\nu}-g_{\mu\nu}(\nabla\Phi)^2/2]/\Phi^2+(\nabla_\mu\nabla_\nu\Phi-g_{\mu\nu}\Box\Phi)/\Phi$ and $\Box\Phi=-4\Lambda_v/(2\omega+3)$.
The chart printer compares every component at points with different $p$, $t_0$, $c$, radius and time, and checks the conformal chart against the cosmic metric pulled back along its time map.
Those checks use no symbolic determinant or simplify operation.
For positive $\rho_v$ and $\varphi_0$, a real nonzero $\chi$ requires $\omega<-3/2$ or $\omega>-5/6$.
Positive Einstein scalar kinetic energy requires $\omega>-3/2$; expansion on the selected $\chi>0$ branch requires $\omega>-1/2$ and acceleration requires $\omega>1/2$.
The selected metric has no observational upper bound on its parameter.
The historical bounds near 25 and 500 constrain phenomenological versions of the model, not the existence of this geometry.
Flat FRW curvature yields $R=6p(2p-1)/(c^2t^2)$, $K=12p^2[(p-1)^2+p^2]/(c^4t^4)$ and vanishing Weyl tensor.

Green and Liddle, Physical Review D 54 (1996), equations (2)-(5), give $g_E=\Phi g_J$ and $\psi=M\sqrt{\omega+3/2}\ln\Phi$ for reduced Planck mass $M$.
The Einstein potential is $V_E=\rho_v\Phi^{-2}=\rho_v\exp[-2\psi/(M\sqrt{\omega+3/2})]$ and $a_E\propto t_E^{(2\omega+3)/4}$.
This conformal change of metric is not an additional coordinate chart of the Jordan geometry.
The canonical exponential-scalar universes retain their independent source and geometry.

For the de Sitter limit, fix $G_0$, $\rho_v$ and elapsed time $\Delta t$ from the reference surface.
Then $t_0^2=c^2(6\omega+5)(2\omega+3)/(32\pi G_0\rho_v)$, $a=(1+\Delta t/t_0)^{\omega+1/2}$ and $\Phi=(1+\Delta t/t_0)^2$.
As $\omega\to\infty$, $\Phi\to1$ and $a\to\exp[\sqrt{8\pi G_0\rho_v/(3c^2)}\,\Delta t]$.
Holding $t$ fixed relative to the singularity while taking that limit would change the specified reference data.

## Drawings

All drawings take $\omega=3/2$, hence $p=2$, and $c=t_0=1$.
This representative accelerating branch is not a claim of observational viability.
Radial null rays have $r\pm c\eta$ constant, with $\eta=-t_0^2/t$ at $p=2$.
Each chart has its own radial spacetime diagram.
The slice $t=t_0$, $\theta=\pi/2$ is a Euclidean plane, embedded as $z=0$ with proper radius $r$.
One moment is drawn, so no movie is needed.
The compactification $T\pm X=2\arctan[(c\eta\pm r)/(ct_0)]$ maps the radial region to $0\le X\le\pi$, $X-\pi\le T<0$.
The regular centre is $X=0$ for $-\pi<T<0$; future infinity is the spacelike boundary $T=0$.
The open sloping past edge $T=X-\pi$ corresponds to $\eta\to-\infty$, $r\to\infty$ with $c\eta+r$ fixed, hence $t\to0$ and divergent curvature.
It is a null singular boundary; the past endpoint of every finite comoving radius is its vertex $(X,T)=(0,-\pi)$.
The spatial-infinity corner $(\pi,0)$ has path-dependent curvature limits and is not claimed to be a regular event.
The observer at the centre has the formal event horizon $r=-c\eta$, of proper radius $ct/(p-1)$, if the background continues indefinitely.
The curvature figure reads the stored $R$ and $K$ and plots their time dependence without recomputing tensors.

## Limits of the approximation

Green and Liddle's equations (11), (13)-(14) and sections III-IV state the constant nucleation approximation and its restrictions.
With $H=p/t$, constant $\Gamma$ implies $\Gamma/H^4\propto t^4$.
The homogeneous solution describes false vacuum regions only.
No bubble wall, collision, reheating or post-transition geometry is inferred from it.
Any future event horizon drawn assumes the homogeneous false vacuum background is continued indefinitely.
