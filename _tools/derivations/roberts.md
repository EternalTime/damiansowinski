# The Roberts solution and critical collapse

A massless scalar field falling in from infinity onto flat space, the same at every scale, with one number $p$ that decides whether it disperses, makes a black hole, or leaves a null singularity.
This note records the five charts, where each comes from, the field equations each is held to, and what each diagram draws.
The sources set $G = c = 1$; the metric file keeps $c$ where a chart has a time, and its null coordinates are lengths.

## Step 1: the double null chart

Oshiro, Nakamura and Tomimatsu's (1), (3), (4) and (5): $ds^2 = -du\,dv + r^2d\Omega^2$ with $r^2 = \tfrac{1}{4}\left((1 - p^2)v^2 - 2uv + u^2\right)$ and $\psi = \pm\tfrac{1}{2}\ln\left(((1 - p)v - u)/((1 + p)v - u)\right)$, for $G_{\mu\nu} = 8\pi T_{\mu\nu}$ with $T_{\mu\nu} = (\partial_\mu\psi\partial_\nu\psi - \tfrac{1}{2}g_{\mu\nu}(\partial\psi)^2)/4\pi$, which is $R_{\mu\nu} = 2\,\partial_\mu\psi\,\partial_\nu\psi$.
Brady's (9) to (11) are the same solution with $ds^2 = -2\,du\,dv$ and $r^2 = \alpha v^2 + \beta u^2 - uv$: at $\beta = 1$, his $u$ is half of this one and $\alpha = (1 - p^2)/4$, so his critical $\alpha = 0$ is $p = 1$.
Burko's (13) and (14) are Oshiro, Nakamura and Tomimatsu's form again, with $p = 2\sigma$ for Roberts's $\sigma$, and his note is the reason the metric file takes this form and not the double null form printed in Roberts's paper of 1989, which solves the field equations only at $\sigma = 0$ and $\sigma = -1$.
The metric file writes the areal radius $R$, since Roberts's own chart has a coordinate $r$ that is not it.
The field is in the region $v > 0$, with flat space before the ray $v = 0$, on which $R = -u/2$ and $\varphi = 0$.
The mass inside a sphere is $-p^2uv/8R$, their (8), which vanishes on $v = 0$ and on $u = 0$.
The Ricci scalar is $8p^2uv/\left(((1 - p)v - u)^2((1 + p)v - u)^2\right)$, which is $p^2uv/2R^4$, their (6), and the Kretschmann scalar is $192p^4u^2v^2$ over the fourth powers of the same two factors, three times the square of the Ricci scalar.

- $p < 1$: $R > 0$ for all $u < 0$, the mass is back to zero on $u = 0$ and nothing crosses that ray, so the spacetime is flat after it, their Figure 2a and Brady's Figure 1.
  On $u = 0$, $R = \sqrt{1 - p^2}\,v/2$, so the flat space after it has the null coordinates $u/\sqrt{1 - p^2}$ and $\sqrt{1 - p^2}\,v$ and its centre is $u = (1 - p^2)v$.
- $p = 1$: $R^2 = u(u - 2v)/4$ vanishes on $u = 0$, a null singularity, where the Kretschmann scalar diverges as $12/u^2v^2$, Brady's case (ii) and Figure 2.
- $p > 1$: the singularity is $u = (1 - p)v$, spacelike, and the apparent horizon, where $\partial_vR = 0$, is $u = (1 - p^2)v$, their (7), of radius $p\sqrt{p^2 - 1}\,v/2$ and mass $p\sqrt{p^2 - 1}\,v/4$.

The published domain is $v > 0$ with $u < 0$ for $p \le 1$ and $u < (1 - p)v$ for $p > 1$.

## Step 2: Roberts's chart

Roberts's (4.5) of 1996, which is (28) of his paper of 1989: $ds^2 = -(1 + 2\sigma)dv^2 + 2\,dv\,dr + r(r - 2\sigma v)d\Omega^2$ with $\varphi = \tfrac{1}{2}\ln(1 - 2\sigma v/r)$.
Burko's (12), $u = (1 + 2\sigma)v - 2r$, carries it to the double null chart, so $r = ((1 + p)v - u)/2$ and $R^2 = r(r - pv)$.
The singularity is $r = pv$, the apparent horizon $r = p(1 + p)v/2$, and the last ray $u = 0$ is $r = (1 + p)v/2$.

## Step 3: the areal radius

Roberts's (4.3) with $\alpha = \sigma v$: his luminosity distance $R^2 = r(r - 2\alpha)$ as the coordinate, with $\lambda^2 = \alpha^2 + R^2$, so $r = \lambda + \alpha$.
Then $dr = \alpha'dv + (\alpha\alpha'dv + R\,dR)/\lambda$ and $ds^2 = -(1 - 2\alpha\alpha'/\lambda)dv^2 + (2R/\lambda)\,dv\,dR + R^2d\Omega^2$, which with $\alpha = pv/2$ is $g_{vv} = -(1 - p^2v/2\lambda)$.
The field is $\tfrac{1}{2}\ln\left((\lambda - \alpha)/(\lambda + \alpha)\right)$, and $u = v - 2\lambda$.
$g^{RR} = \lambda(2\lambda - p^2v)/2R^2$ vanishes on the apparent horizon, $R = p\sqrt{p^2 - 1}\,v/2$.
The chart names $\lambda$, and `print_charts.roberts_lambda` writes every value in $\lambda$, $v$ and $p$ with at most one $R$ in front.

## Step 4: the time and radius chart

Roberts's (4.6) and (4.7): with $v' = \sqrt{1 + 2\sigma}\,v$, $r' = r/\sqrt{1 + 2\sigma}$ and $t' = v' - r'$, $ds^2 = -dt'^2 + dr'^2 + r'(r' - 2\sigma t')d\Omega^2$.
The metric file writes $t$ for $t'$ and $\rho$ for $r'$, so $u = \sqrt{1 + p}\,(ct - \rho)$ and $v = (ct + \rho)/\sqrt{1 + p}$.
The plane of $t$ and $\rho$ is flat, and the whole collapse is in $R = \sqrt{\rho(\rho - pct)}$.
The first ray is $ct = -\rho$, the last ray $ct = \rho$, the singularity $ct = \rho/p$, and the apparent horizon $ct = (2 - p)\rho/p$, from $(2\rho - pct)^2 = p^2\rho^2$.
The chart is called `diagonal`.

## Step 5: the scaling chart

Frolov's (10) to (12) for his $dS^2 = -2\,du\,dv$: $x = \tfrac{1}{2}\ln(1 - v/u)$ and $s = -\ln(-u)$, in which the critical solution is $R = e^{x - s}$ and $\Phi = x$.
With a length restored and the null coordinates of Step 1, $u = -2\ell e^{-\tau}$ and $v = \ell e^{-\tau}(e^{2x} - 1)$, his $s$ written $\tau$ since a coordinate $s$ would make the line element's $ds^2$ ambiguous.
Then $ds^2 = 2\ell^2e^{2(x - \tau)}\left((1 - e^{-2x})d\tau^2 - 2\,d\tau\,dx\right) + R^2d\Omega^2$ with $R^2 = \ell^2e^{2(x - \tau)}(\cosh^2x - p^2\sinh^2x)$, which is Frolov's (12) and (13) at $p = 1$ and $\ell = 1$.
The field is $\tfrac{1}{2}\ln\left((1 - p\tanh x)/(1 + p\tanh x)\right)$, a function of $x$ alone, and the only component of the Ricci tensor is $R_{xx} = 2p^2/(\cosh^2x - p^2\sinh^2x)^2$.
The singularity is $x = \tanh^{-1}(1/p)$ and the apparent horizon $x = \tanh^{-1}(1/p^2)$, each a line of constant $x$.
`print_charts.roberts_scaling` writes every value in $e^{2x}$ and $e^{2\tau}$ with the two factors of the areal radius paired as $4e^{2x}(\cosh^2x - p^2\sinh^2x)$.

## Step 6: what is checked

`print_charts.roberts_check` holds every chart, at six random points of its domain in forty digits, to $R_{\mu\nu} = 2\,\partial_\mu\varphi\,\partial_\nu\varphi$ in every slot, to the wave equation, to flat space at $p = 0$, and to being the double null chart pulled back.
Not published as a chart: Oshiro, Nakamura and Tomimatsu's diagonal chart of $t$ and the areal radius, their (11) to (14), whose time is given only as an integral.

## Step 7: the diagrams

Every diagram is drawn for $p = 9/10$, $1$ and $2$, three spacetimes of one line element.
Nothing in the solution sets a scale, so every length is in units of any length $\ell$.
The dispersing case is drawn at $p = 9/10$ because a weaker field bends a moment of space too little to show: at $p = 1/2$ the embedded surfaces rise by under a twentieth of their width.

The spacetime diagrams draw the double null chart against $(v - u)/2$ and $(u + v)/2$, Roberts's chart against $r$ and $(1 + p)v - r$, which are $\sqrt{1 + p}$ times $\rho$ and $ct$, the areal chart against $R$ and $(1 + p/2)v - R$, and the scaling chart against $x$ and $\tau + x$, a time function since $g^{\tau\tau} = 0$ and $g^{\tau x}$, $g^{xx}$ are negative.
Where the drawn radius is a null coordinate, as $v$ is in the double null chart and $x$ at fixed $\tau$ in the scaling chart, the curve on which the areal radius is stationary along it is the apparent horizon, which the row option `null_radius` marks.
The areal chart degenerates on $R = 0$, where $g_{vR} = R/\lambda$ vanishes, so its rays at $p = 2$ stop at the edge of the domain.
`null_rays.py --verify` checks the rays of every view against $u$ and $v$.

The conformal diagrams bring the null coordinates in by $\arctan(u/\ell)$ and $\arctan(v/\ell)$.
For $p < 1$ the centre of the flat space after the last ray is $u = (1 - p^2)v$, so $u$ is divided by $1 - p^2$ there, which is continuous across $u = 0$ and puts that centre on $X = 0$: the whole is Minkowski's triangle.
At $p = 1$ the singularity is the null line from the centre to $i^+$, and for $p > 1$ a spacelike line from the centre to $i^0$, level at $p = 2$ where $u = -v$.

The embedding diagrams are the equator of a moment of Roberts's $t$, read in the diagonal chart, where the metric on it is $d\rho^2 + \rho(\rho - pct)d\phi^2$.
$(dR/d\rho)^2 - 1 = p^2c^2t^2/4R^2$, so no moment but $t = 0$ stands in flat space and every one stands in three dimensional Minkowski space, with $dZ/d\rho = pc|t|/2R$ and $Z = pct\ln\left(\sqrt{\rho} + \sqrt{\rho - pct}\right)$ up to a constant.
Each moment is drawn with its rim $\rho = 5\,\ell$ at $Z = 0$, from the first ray $\rho = -ct$ before $t = 0$, and after it from the last ray $\rho = ct$ for $p < 1$ and from the singularity $\rho = pct$ for $p \ge 1$.
Toward the singularity $dZ/dR = pct/(2\rho - pct) \to 1$ and $1 - (dZ/dR)^2 = 4R^2/p^2c^2t^2$, so the profile starts $2 \times 10^{-4}$ of $pct$ short of it, steps by a fiftieth of $R$, and is written to fourteen decimals.
The moments marked are $ct = \pm 0.8\,\ell$ and $\pm 1.6\,\ell$, and the movie runs from $-1.6\,\ell$ to $1.6\,\ell$ with a frame every $\ell/10$.
