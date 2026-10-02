# Witten's black hole in two dimensions

The black hole of dilaton gravity in two dimensions, the solution of the string's low energy equations that Witten found as an exact conformal field theory and that Mandal, Sengupta and Wadia found from the equations themselves.
This note records the seven charts, where each comes from, the field equations each is held to, and what each diagram draws.
The sources write $G = c = \hbar = 1$; the metric file keeps $c$.
Its seven charts are written by `_tools/derivations/print_charts.py --metric witten_black_hole`, and `verify_metrics.py --system witten_black_hole/<chart>` checks each in a second.

Equation numbers are those of the preprints: Witten's IASSNS-HEP-91/12, the text of Physical Review D 44, 314; Mandal, Sengupta and Wadia's IASSNS-HEP-91/10; Callan, Giddings, Harvey and Strominger's arXiv:hep-th/9111056; Horowitz's arXiv:hep-th/9210119; Strominger's arXiv:hep-th/9501071; and Grumiller, Kummer and Vassilevich's arXiv:hep-th/0204253.

## Step 1. Two dimensions

It is the first spacetime of two dimensions in the collection, and two things about its tensors follow from the dimension alone.

The Riemann tensor of any metric in two dimensions is $R_{\mu\nu\rho\sigma} = \tfrac{1}{2}R\,(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho})$, all trace.
So the Einstein tensor vanishes identically, the Kretschmann scalar is $K = R^2$, and the Weyl tensor, what is left of the Riemann tensor with its traces removed, is zero.
The formula the checker removes the traces by divides by $n - 2$, so `Geometry.weyl_llll` in `verify_metrics.py` returns zero for $n < 3$ and never evaluates it.
Every chart publishes an Einstein tensor and a Weyl tensor with no component, as the BTZ black hole publishes its Weyl tensor, and `witten_black_hole_check` holds $G_{\mu\nu} = 0$ and $K = R^2$ in each.

The Einstein equation is therefore no equation here, and the metric is fixed by the dilaton's.

## Step 2. The field equations

The action is Callan, Giddings, Harvey and Strominger's (1) without its matter, $S = \frac{1}{2\pi}\int d^2x\sqrt{-g}\,e^{-2\Phi}\left(R + 4(\nabla\Phi)^2 + 4\lambda^2\right)$, with $\lambda$ an inverse length.
Its field equations are Strominger's (3.7) and (3.8), which together are

$$R_{\mu\nu} + 2\nabla_\mu\nabla_\nu\Phi = 0,\qquad R + 4\lambda^2 + 4\nabla^2\Phi - 4(\nabla\Phi)^2 = 0.$$

`witten_black_hole_check` in `print_charts.py` holds every chart to both, slot by slot, with the chart's own dilaton, before the chart is written.
Witten's dilaton has the other sign and twice the size, his (12), $\Phi_{\text{W}} = 2\ln\cosh r + \text{const}$; the collection takes the sign and size of the action above, in which $e^{\Phi}$ is the string coupling.

The empty solution is flat space with the dilaton linear, $\Phi = -\lambda x$, the linear dilaton vacuum.

## Step 3. The charts

`witten`, his (22) with $\lambda$ written out: $ds^2 = -\tanh^2(\lambda r)\,c^2dt^2 + dr^2$ with $\Phi = \Phi_0 - \ln\cosh(\lambda r)$.
$r$ is the proper distance from the horizon $r = 0$, and no mass appears in the metric: by his (48) and (49) the mass is set by the additive constant of the dilaton, the value $\Phi_0$ on the horizon.
$R = 4\lambda^2/\cosh^2(\lambda r)$, which is $4\lambda^2$ on the horizon and falls off exponentially.

`schwarzschild_gauge`, Mandal, Sengupta and Wadia's first method: the gauge in which the dilaton is proportional to a coordinate, $ds^2 = -f\,c^2dt^2 + dx^2/f$ with $f = 1 - m\,e^{-2\lambda x}$ and $\Phi = -\lambda x$.
Theirs is $g = 1 - a\,e^{Q\eta}$ with $\Phi = Q\eta/2$; their note (i) allows either sign of $Q$, and $Q = -2\lambda$, $\eta = x$, $a = m$ puts the flat end at $x \to \infty$ and the singularity at $x \to -\infty$.
The parameter $m = e^{-2\Phi_0}$ is dimensionless and positive, and it is the $M/\lambda$ of Callan, Giddings, Harvey and Strominger, whose (11) has $e^{-2\phi} = M/\lambda$ on the horizon.
A shift of $x$ changes $m$, so one value of $m$ is every other in another origin of $x$; the horizon is $x = \ln(m)/2\lambda$.
The map from Witten's chart is $e^{2\lambda x} = m\cosh^2(\lambda r)$, so $f = \tanh^2(\lambda r)$ and $dx = \tanh(\lambda r)\,dr$.
$R = 4\lambda^2m\,e^{-2\lambda x}$.

`dilaton`, Horowitz's (4.20): $e^{-2\Phi}$ itself as the coordinate, $w = e^{-2\Phi} = e^{2\lambda x}$, in which $ds^2 = -(1 - m/w)\,c^2dt^2 + dw^2/4\lambda^2w(w - m)$.
His $r$ is $w$ up to a constant, his $M$ is $m$ in that unit, and his $k$ is $1/\lambda^2$.
The horizon is $w = m$ and the singularity $w = 0$, where $R = 4\lambda^2m/w$ diverges and the string coupling $e^{\Phi}$ with it.

`conformal`, Witten's (25) and the others' (15): the tortoise coordinate $\sigma$, $d\sigma = dx/f$, so $e^{2\lambda\sigma} = e^{2\lambda x} - m = m\sinh^2(\lambda r)$ and $ds^2 = (-c^2dt^2 + d\sigma^2)/(1 + m\,e^{-2\lambda\sigma})$.
It covers the outside of the horizon, which is $\sigma \to -\infty$.

`kruskal`, Witten's (26) to (28): $V = \sinh(\lambda r)\,e^{\lambda ct}$ and $U = -\sinh(\lambda r)\,e^{-\lambda ct}$, dimensionless, with $ds^2 = -dU\,dV/\lambda^2(1 - UV)$.
$UV = 1 - e^{-2\Phi}/m$ in every chart: $-\sinh^2(\lambda r)$, $1 - e^{2\lambda x}/m$, $1 - w/m$ and $-e^{2\lambda\sigma}/m$.
Mandal, Sengupta and Wadia's second method gives the same form, and Callan, Giddings, Harvey and Strominger's (11) is it with $x^\pm = \sqrt{m}\,(V, U)/\lambda$.
The chart names $V$ the advanced coordinate and $U$ the retarded one, as Kruskal's are named for Schwarzschild, so the outside of the hole is $U < 0 < V$.
The horizons are $UV = 0$ and the singularity is $UV = 1$, where $R = 4\lambda^2/(1 - UV)$ diverges.
Witten's regions V and VI, beyond $UV = 1$, where the metric changes sign and "time flows sideways", are no part of the spacetime published: the domain is $UV < 1$.

`eddington_finkelstein_ingoing` and `eddington_finkelstein_outgoing`, the gauge of Grumiller, Kummer and Vassilevich's (3.26) with the Killing norm of their (3.74): $ds^2 = -f\,dv^2 + 2\,dv\,dx$ and $ds^2 = -f\,du^2 - 2\,du\,dx$, with $v, u = ct \pm \sigma$ lengths and $f$ and $\Phi$ those of the Schwarzschild gauge.
Their asymptotic region is at $r \to -\infty$, the other sign of the coordinate.

Each chart after the first is checked to be the first: Witten's metric carried along the map into the chart is the chart's published metric, slot by slot.

## Step 4. The temperature and the time inside

The surface gravity is $\kappa = c^2f'(x_h)/2 = \lambda c^2$ for every $m$, so the Hawking temperature is $\hbar c\lambda/2\pi k_B$ whatever the mass, Strominger's $T = \lambda/2\pi$.
In Witten's chart the Euclidean section is $dr^2 + \tanh^2(\lambda r)\,d\theta^2/\lambda^2$ with $\theta = i\lambda ct$, his (9), and it is smooth at $r = 0$ only for the period $2\pi$ of $\theta$, a period $2\pi/\lambda c$ of the Euclidean time.

Inside the horizon $x$ is the time, and the proper time along a curve of constant $t$ from the horizon to the singularity is $\int dx/c\sqrt{m\,e^{-2\lambda x} - 1} = \pi/2\lambda c$, with $y = \sqrt{m}\,e^{-\lambda x}$ the integral of $dy/\lambda cy\sqrt{y^2 - 1}$ from $1$ to $\infty$.
Those curves are geodesics, and a timelike curve between two surfaces of constant $x$ is longest when it stays at constant $t$, since $d\tau^2 = dx^2/|f| - |f|\,c^2dt^2$ there; so no observer spends longer than $\pi/2\lambda c$ between the horizon and the singularity, as the caption of the Schwarzschild gauge's plane says.
`_tools/test_witten_black_hole.py` holds both numbers.

## Step 5. The diagrams

Every diagram is drawn at $\lambda = 1$, the unit of length being $1/\lambda$, and $m = 1$, where the horizon is $x = 0$ and $w = 1$.

The spacetime diagrams draw the plane of each chart, seven views.
The tortoise coordinate is $\ln\sinh r$, $\tfrac{1}{2}\ln|e^{2x} - 1|$, $\tfrac{1}{2}\ln|w - 1|$ and $\sigma$ itself, and `null_rays.py --verify` checks the rays of each view against $ct \pm$ it, and those of the Kruskal plane against $U$ and $V$.
In the dilaton chart $dw/d(ct) = \pm 2\lambda(w - m)$, so the rays are exponentials and meet the singularity $w = 0$ at a finite slope.
The Kruskal chart's domain is the inequality $UV < 1$, which no interval of one coordinate states, so its row declares `where="1 - U*V"`: the view is hatched where that is not positive, and no ray or cone is drawn there.
The singularity is declared in `singular_zero` and checked on the Kretschmann scalar.

The conformal diagram is Kruskal's hexagon, $p = \arctan U$ and $q = \arctan V$, with the singularity $UV = 1$ on the straight lines $T = \pm\pi/2$ since $\tan(p + q) = (U + V)/(1 - UV)$.
It has one view for each chart, each tinting what its chart covers, and `conformal.py` checks every map against the published metric of its chart, and one event to land on one point through all seven.
Each point of the diagram is a single event.

The embedding diagram is the cigar, the Euclidean section.
A moment of this spacetime is a line, which has nothing to embed, and the surface the literature draws is Witten's figure 1.
`Slice` reaches it with a `swept` that turns the time through the imaginary direction, $t = -i\phi$: $g_{\phi\phi}$ is then $-g_{tt} = \tanh^2 r$ and the slice is a surface of revolution with $\rho = \tanh r$ and $d\rho/ds = 1/\cosh^2 r \le 1$, so it embeds in flat space everywhere.
Its height is $z = \operatorname{arsinh}(\cosh r) - \sqrt{1 + 1/\cosh^2 r} + \sqrt{2} - \operatorname{arsinh} 1$, from $dz/dr = \tanh r\sqrt{1 + 1/\cosh^2 r}$, which `embedding.py` checks.
The Euclidean and Lorentzian sections meet along the moment $t = 0$: under $t = -i\theta/\lambda c$ the Kruskal coordinates are $V = \sinh(\lambda r)\,e^{-i\theta}$ and $U = -\sinh(\lambda r)\,e^{i\theta}$, real at $\theta = 0$, the right exterior's $t = 0$, and at $\theta = \pi$, the left exterior's.
So the cigar's two meridians $\theta = 0$ and $\pi$ are the line $U + V = 0$ through the bifurcation point, and that is the moment `slices.py` marks on the other diagrams: outside the horizon in each chart, and on both sides of it on the Kruskal plane and the conformal views.
The geometry is static, so it is one surface and no movie.
