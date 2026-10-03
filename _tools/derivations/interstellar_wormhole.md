# Interstellar's wormhole

The Dneg wormhole of Oliver James, Eugénie von Tunzelmann, Paul Franklin and Kip S. Thorne, "Visualizing Interstellar's Wormhole", Am. J. Phys. 83, 486 (2015), doi:10.1119/1.4916949, arXiv:1502.03809, read in full from the arXiv source (v3).
Every equation number below is theirs.

## Step 1: the metric

Their (1) is the general wormhole metric, $ds^2 = -dt^2 + d\ell^2 + r(\ell)^2(d\theta^2 + \sin^2\theta\,d\phi^2)$ in units $G = c = 1$, with $\ell$ the proper distance along a radial line.
The Dneg wormhole without gravity is (1) with (5):

- (5c): $r = \rho$ for $|\ell| < a$, the cylinder of length $2a$;
- (5a): $r = \rho + \mathcal{M}[x\arctan x - \tfrac{1}{2}\ln(1 + x^2)]$ for $|\ell| > a$;
- (5b): $x = 2(|\ell| - a)/\pi\mathcal{M}$.

The page writes their $\mathcal{M}$ as $M$, since the checker's reader takes no calligraphic letters, and restores $c$: $-c^2dt^2$.
Their (8) adds a potential, $-(1 + 2\Phi)dt^2$ with $|\Phi| \lesssim 10^{-12}$, which they drop for every ray they trace; the page publishes the metric without it.

## Step 2: the three charts

All three are their coordinates $(t, \ell, \theta, \phi)$.

- `proper_distance` is (1) over the whole line of $\ell$ with $r(\ell)$ a declared function, so its values hold for every $r(\ell)$ and for (5) in particular, as Morris and Thorne's proper radial chart does with $\Phi = 0$.
- `cylinder` is $|\ell| \le a$ with $r = \rho$: the product of the flat plane of $t$ and $\ell$ with a sphere of radius $\rho$, which is Plebański and Hacyan's flat plane times a sphere with $b = \rho$.
- `flare` is $\ell \ge a$ with (5a) written out through the defined names $x$ and $r$; the flare $\ell \le -a$ is its mirror image $\ell \to -\ell$.

No chart covers both flares with the explicit radius: written with $|\ell|$, the checker's derivative of $\mathrm{sgn}(\ell)$ puts a $\delta(\ell)$ at $\ell = 0$, inside the cylinder, where the flare's formula does not hold.

## Step 3: the curvature

For $-dt^2 + d\ell^2 + r^2d\Omega^2$ the frame components of the Riemann tensor are $-r''/r$ on the planes holding $\ell$ and $(1 - r'^2)/r^2$ on the sphere, so
$K = 8r''^2/r^2 + 4(1 - r'^2)^2/r^4$ and $R = 2(1 - r'^2)/r^2 - 4r''/r$.
On the flare, their footnote 19 gives $r' = \tfrac{2}{\pi}\arctan x$, and then $r'' = 4/(\pi^2M(1 + x^2))$.
`interstellar_wormhole_check` in `print_charts.py` holds every chart to these, and the flare to $r = \rho$ and $r' = 0$ at $\ell = a$ and $r' \to 1$ far away.
$r$ and $r'$ are continuous at the mouths and $r''$ jumps there from $0$ to $4/\pi^2M$, so the curvature jumps across each mouth and carries no delta.
On the cylinder $K = 4/\rho^4$.

The flare's values are printed in $x$, $\arctan(x)$, $r$, $1 + x^2$ and $\pi^2 - 4\arctan(x)^2 = \pi^2(1 - r'^2)$.
`norm` in `verify_metrics.py` factors the argument of every arctangent, so that $\arctan(2(\ell - a)/\pi M)$ and $\arctan(2a/\pi M - 2\ell/\pi M)$ are one generator and its negative; before that change the two spellings were two generators and no value of the flare read back.

## Step 4: the lensing width

Their (6) embeds the equator: $dz/d\ell = \sqrt{1 - r'^2}$.
The surface stands at $45°$ where $r' = 1/\sqrt{2}$, at $x = \tan(\pi/2\sqrt{2})$, and there $r - \rho = W$ with
$W/M = x\arctan x - \tfrac{1}{2}\ln(1 + x^2) = (\pi/2\sqrt{2})\tan(\pi/2\sqrt{2}) - \ln\sec(\pi/2\sqrt{2}) = 1.42953$,
their (7).
Their footnote 19 prints $1.42053$, a slip: its own closed form gives $1.42953$, and their figure 3 has $W/\rho = 0.715$ at $M/\rho = 0.5$.
`embedding.py` checks the number on the drawn surface.

## Step 5: the drawings

Every diagram is drawn at $a = \rho$ and $M = \rho/2$, the wormhole of their figure 3, in units of $\rho$; `null_rays.iw_radius` is (5) with three derivatives, declared to sympy as `iw_r`.
The plane of $t$ and $\ell$ is $-c^2dt^2 + d\ell^2$ in every chart, so the rays keep $ct \mp \ell$ and the proper distance is its own tortoise coordinate: the conformal diagram is the full diamond, as for Ellis's wormhole, with the cylinder the band between the mouths.
The film's own wormhole, $2a = 0.01\rho$ and $W = 0.05\rho$, is the History's; at that size the flare is a fiftieth of the radius, and the figure-3 wormhole shows the cylinder and the flares at a scale a drawing can carry.
