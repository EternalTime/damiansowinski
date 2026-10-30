# Moncrief generalised Taub-NUT cosmologies

## Scope

The subject is the infinite-dimensional family of toroidal vacuum cosmologies introduced in 1981.
The primary bibliographic records are verified in the task literature ledger.
The representation retains both Gowdy wave amplitudes and restricts them to analytic solutions with a compact horizon.
The original full texts remain unread; no body quotation or equation is attributed to them.
The equation-bearing sources are Chruściel, Isenberg, and Moncrief (1990), equations (1)-(4) and section 4, and Chruściel and Lake (2004), equations (1.4), (1.5), (3.7), and (4.2).

## Defining equations

Use dimensionless coordinates $(t,\theta,\psi,\phi)$, a length $L>0$, and a positive constant $\kappa$.
The three spatial coordinates have period $2\pi$.
Set $P=\ln t+2b$, $Q=q$, and $\gamma=-\ln t-4b-4f-4\ln\kappa$ in Chruściel and Lake's metric (1.4), reversing their negative areal time.
The resulting metric is
$$ds^2=L^2\left\{\kappa^2e^{2b+2f}(-dt^2+d\theta^2)+t^2e^{2b}(d\psi+q\,d\phi)^2+e^{-2b}d\phi^2\right\}.$$
Their equations (4.2) and (1.5) become
$$b_{tt}+t^{-1}b_t-b_{\theta\theta}=\tfrac12t^2e^{4b}(q_t^2-q_\theta^2),$$
$$q_{tt}+3t^{-1}q_t-q_{\theta\theta}+4(b_tq_t-b_\theta q_\theta)=0,$$
$$f_t=t(b_t^2+b_\theta^2)+\tfrac14t^3e^{4b}(q_t^2+q_\theta^2),\qquad f_\theta=2tb_tb_\theta+\tfrac12t^3e^{4b}q_tq_\theta.$$
These are a change of variables in the sourced equations, not an independent tensor computation.
The page must specify analytic functions of $t^2$ and $\theta$ near $t=0$, periodicity, and $f(0,\theta)=0$.
The definition includes only solutions satisfying these equations and conditions.
No claim that arbitrary smooth horizon data admit a vacuum continuation follows from this definition.
The analytic polarised sector is infinite dimensional by the freely chosen analytic function in the 1990 paper's section 4.
Its branches $\pi=0$ and $\pi=1$ are exchanged by interchanging the two Killing coordinates.

## Horizon extension

Let $T=t^2$ and $\chi_s=\psi+s\kappa\ln t$, with $s=+1$ or $s=-1$.
Write the analytic continuations of $b,f,q$ as $B(T,\theta),F(T,\theta),Q(T,\theta)$.
Define $H=\kappa^2(e^{2F}-1)/(4T)$ by its analytic value at $T=0$.
Direct substitution gives
$$ds^2=L^2\left\{e^{2B}\left[T(d\chi_s+Qd\phi)^2-s\kappa\,dT(d\chi_s+Qd\phi)-H\,dT^2\right]+\kappa^2e^{2B+2F}d\theta^2+e^{-2B}d\phi^2\right\}.$$
This is a derived chart, to be checked by a numerical pullback in the chart specification.
Both signs must be retained as distinct extension choices, with their compact identifications stated.
The determinant at the horizon is $-L^8\kappa^4e^{4B(0,\theta)}/4$.
The periodic generator $\partial_{\chi_s}$ has squared norm $L^2Te^{2B}$ and is timelike for $T<0$.
Analyticity supplies a vacuum continuation only within the continuation's domain of convergence.
The extension domain must be local unless a particular member's functions are continued farther.

## Drawings

Use $q=0$, $\kappa=1$, $L=1$, and $b=\epsilon J_0(t)\cos\theta$, with $\epsilon=1/4$.
The corresponding constraint solution is
$$f=\epsilon^2\left[\tfrac12t^2(J_0(t)^2+J_1(t)^2)-tJ_0(t)J_1(t)\cos^2\theta\right].$$
The Bessel identities make both first-derivative constraints hold, and $f(0,\theta)=0$.
The Bessel functions are even analytic functions of $t$ in these combinations, so the extension is defined by power series in $T$, including negative $T$.
This member is an explicitly derived inhomogeneous sector of the field equations.
It is not the full family's parameterization.
A spacetime drawing should include the compact collapsing direction and the regular horizon chart; a plot only in $(t,\theta)$ would hide the extension mechanism.
An embedding may use a cut surface with identified ends if its induced metric is checked and the identification is labelled.
The embedding must not claim to embed the entire compact torus as that cut surface.
The conformal diagram must specify the two-dimensional quotient or section and its periodic identification.
The curvature figure must use the stored scalar, with no new frame kind.

## Checks

The chart specifications and dedicated spacetime check remain to be implemented after History approval.
The pointwise checks must exercise nonzero $q$ derivatives as well as the polarised drawing member.
The constraints, regular horizon determinant, sign change of the compact generator, and pullbacks of both extension signs are distinct checks.
The page-scoped checker runs once after the page is complete.

## Implemented charts and figures

All three charts retain both wave amplitudes in their tensor expressions.
The printer checks the vacuum equations at generic jets with nonzero derivatives of the second amplitude, and compares both logarithmic coordinate pullbacks at points.
A named positive radius $\delta$ restricts each extension to its local analytic domain.
The dedicated spacetime test evaluates the horizon determinant numerically after writing out the removable quotient, and checks the sign of the periodic generator on both sides.

The spacetime drawings show the totally geodesic section $\theta=\phi=0$ of the polarised Bessel member, with the edges of each periodic angle identified.
On that section, derivatives of every metric coefficient transverse to the section vanish.
The conformal generator also checks this against the published Christoffel symbols.
For the regular section, the entire functions $A(T)=J_0(\sqrt{T})$ and $D(T)=J_1(\sqrt{T})/\sqrt{T}$ have real power series on either side of the horizon.
The constraint reads $F_T=T D(T)^2/32$ there, so evaluating $F/T$ by its convergent series avoids cancellation in $H$ at zero.
The numerical drawing functions retain forty terms of this convergent series on $|T|\le2$.
The spacetime is still defined by the analytic functions and the vacuum equations, with the finite series used solely to evaluate the drawing coefficients.

For the conformal map define $w(T)=\int_0^T(e^{F(u)}-1)/(2u)\,du$.
The plus extension has $U=e^{\chi+w}$ and $V=T e^{-\chi+w}$, while the minus extension has $U=T e^{\chi+w}$ and $V=e^{-\chi+w}$.
In both cases the induced line element is $-L^2e^{2B-2w}dU\,dV$ at $\kappa=1$.
The drawn coordinates are $p=\arctan U$ and $q=\arctan V$.
The periodic identification is $(U,V)\mapsto(e^{2\pi}U,e^{-2\pi}V)$.
Finite cuts at $T=\pm1$ bound the local views; no global analytic extension or conformal infinity is inferred from them.
The embedded moment is $t=T=1$, where both angular shifts vanish, and its compact circle is marked in every spacetime and conformal view.

The embedding is the torus $t=1$, $\phi=0$, cut along $\theta=0$ and $2\pi$ with its ends identified.
Its radius is $Lte^b$, and its height obeys $z_\theta^2=L^2\kappa^2e^{2b+2f}-L^2t^2e^{2b}b_\theta^2$.
The generator checks that this is positive throughout the selected slice and that the induced metric is preserved.
The open tube with its stated identification represents the cut torus.
A single moment is drawn, so no movie or stack registration is required.
The curvature figure reads only the stored Kretschmann scalar on $0.2\le t\le4$ over one full spatial period.
All four kinds of diagram are present, so no absent-diagram registration is needed.
