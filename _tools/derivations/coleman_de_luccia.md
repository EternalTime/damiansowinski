# The Coleman-De Luccia bubble

The five charts of `coleman_de_luccia.json` are written by `print_charts.py --metric coleman_de_luccia`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records where each chart comes from, what `coleman_de_luccia_check` holds it to, and what the diagrams draw.
Equation numbers are those of Coleman and De Luccia, Phys. Rev. D 21, 3305 (1980), read in the preprint SLAC-PUB-2463, and of Bucher, Goldhaber and Turok, Phys. Rev. D 52, 3314 (1995), read in arXiv:hep-ph/9411206.

## Step 1. The chart outside the light cone

Coleman and De Luccia assume the bounce is invariant under rotations of Euclidean four-space, which leaves the metric $d\xi^2 + \rho(\xi)^2d\Omega_3^2$, their (3.1), with one unknown function.
Continued to real time the rotations become the Lorentz group about the event the bubble is centred on, and outside that event's light cone the metric is their (4.1), $d\xi^2 + \rho(\xi)^2$ times the metric of the unit hyperboloid with spacelike normal, de Sitter space of three dimensions.
Their signature is $(+,-,-,-)$ and they leave the hyperboloid without coordinates; Bucher, Goldhaber and Turok write it out in their (5.4), and with the site's signature and the name $\psi$ for the rapidity it is

$$ds^2 = d\xi^2 + \rho^2\left(-d\psi^2 + \cosh^2\psi\,d\theta^2 + \cosh^2\psi\sin^2\theta\,d\phi^2\right).$$

The chart keeps $\rho$ free, so it holds for a wall of any thickness.
Its published $G^\xi{}_\xi = 3\left((\partial_\xi\rho)^2 - 1\right)/\rho^2$ set equal to $\kappa T^\xi{}_\xi = \kappa\left(\tfrac12\phi'^2 - U\right)$ is their (3.4), $\rho'^2 = 1 + \tfrac{\kappa}{3}\rho^2\left(\tfrac12\phi'^2 - U\right)$.
In a vacuum of cosmological constant $\Lambda$ that gives $\rho'^2 = 1 - \Lambda\rho^2/3$, whose solutions through $\rho(0) = 0$ are $\xi$, $\ell\sin(\xi/\ell)$, their (4.5), and $\ell\sinh(\xi/\ell)$, their (4.11), with $\ell = \sqrt{3/|\Lambda|}$; the check substitutes each and finds $G_{\mu\nu} = -\Lambda g_{\mu\nu}$ in every slot, and the Riemann tensor zero for $\rho = \xi$.
$\psi$ is a pure number and is itself the chart's $x^0$.

## Step 2. The thin wall

A thin wall of surface energy density and tension $\sigma$ on the hyperboloid $\rho = \bar\rho$ has $T^\psi{}_\psi = -\sigma\,\delta(\xi - \xi_w)$.
The published $G^\psi{}_\psi = \left(2\rho\,\partial_\xi^2\rho + (\partial_\xi\rho)^2 - 1\right)/\rho^2$ is singular there only through $\partial_\xi^2\rho$, so Einstein's equation gives Israel's condition in the form

$$\partial_\xi\rho\,\big|_{\text{inside}} - \partial_\xi\rho\,\big|_{\text{outside}} = \frac{4\pi G\sigma\bar\rho}{c^4}.$$

With $\partial_\xi\rho = \sqrt{1 - \Lambda\rho^2/3}$ on each side this is one equation for $\bar\rho$.
Coleman and De Luccia's $\bar\rho_0 = 3\sigma/\epsilon$ and their length $\Lambda = (\kappa\epsilon/3)^{-1/2}$, written $\ell$ here, make $4\pi G\sigma/c^4 = \bar\rho_0/2\ell^2$, and the equation is solved by their (3.15), $\bar\rho = \bar\rho_0/(1 + (\bar\rho_0/2\ell)^2)$, for a decay into flat space and by their (3.18), $\bar\rho = \bar\rho_0/(1 - (\bar\rho_0/2\ell)^2)$, for a decay of flat space, which `coleman_de_luccia_junction` checks.
They found those radii by making the action stationary; the junction condition gives the same radii from the metric alone.

The drawings take $\bar\rho_0 = \ell$ in both cases.
Decay into flat space: $\bar\rho = 4\ell/5$, and $\partial_\xi\rho$ drops from $1$ to $3/5$, so outside the wall $\rho = \ell\sin((\xi - \xi_0)/\ell)$ with $\xi_0 = 4\ell/5 - \ell\arcsin(4/5)$.
Decay of flat space: $\bar\rho = 4\ell/3$ at $\xi = \ell\ln 3$, and $\partial_\xi\rho$ drops from $5/3$ to $1$, so outside the wall $\rho = \xi - \ell\ln 3 + 4\ell/3$.

## Step 3. The open universe

Through the light cone, $\xi = ic\tau$ and $\rho(ic\tau) = i\,a(\tau)$, and the hyperboloid with spacelike normal becomes the one with timelike normal, hyperbolic space, through $\psi = \chi + i\pi/2$.
That is their (4.2), "a Robertson-Walker universe of open type", and Bucher, Goldhaber and Turok's (5.1):

$$ds^2 = -c^2d\tau^2 + a^2\left(d\chi^2 + \sinh^2\chi\,d\theta^2 + \sinh^2\chi\sin^2\theta\,d\phi^2\right).$$

The check makes that continuation on the printed wall chart and compares every slot.
The vacuum scale factors are $c\tau$, $\ell\sinh(c\tau/\ell)$ and $\ell\sin(c\tau/\ell)$, their (4.12) for the last, which returns to zero at $c\tau = \pi\ell$.
The conformal time chart is the same universe with $a\,d\eta = c\,d\tau$, checked by pulling the open chart back; in the vacua $\eta = \ln(c\tau/\ell)$, $\ln\tanh(c\tau/2\ell)$ and $\ln\tan(c\tau/2\ell)$, so $a = \ell e^\eta$, $-\ell/\sinh\eta$ and $\ell/\cosh\eta$.

## Step 4. The static charts

On either side of a thin wall the spacetime is a vacuum, and its static chart about the bubble's centre is the chart Blau, Guendelman and Guth and Berezin, Kuzmin and Tkachev follow a wall in.
With $H^2 = \Lambda/3$ of either sign, the wall chart of that vacuum is the static chart pulled back along

$$r = \rho\cosh\psi, \qquad \tanh(Hct) = \frac{H\rho\sinh\psi}{\partial_\xi\rho},$$

which is $ct = \xi\sinh\psi$ at $\Lambda = 0$ and holds a tangent for $\Lambda < 0$; the check pulls back at random points in forty digits for each sign.
On the wall, $\rho = \bar\rho$, eliminating $\psi$ gives

$$r_w^2 = \bar\rho^2 + \left(1 - \frac{\Lambda\bar\rho^2}{3}\right)\frac{\tanh^2(Hct)}{H^2},$$

the hyperbola $\bar\rho^2 + c^2t^2$ in flat space.
The check holds that curve to carrying the metric $\bar\rho^2(-d\psi^2 + \cosh^2\psi\,d\Omega^2)$.
In de Sitter space the wall reaches the horizon $r = \ell$ only as $t \to \infty$, and in anti-de Sitter space it reaches infinite $r$ at $ct = \pi\ell/2$.

## Step 5. What the diagrams draw

The spacetime diagrams draw the two bubbles of Step 2 in the wall chart and in the static charts on both sides, the open universe for each sign of the new vacuum's energy, and the conformal time chart for the negative one.
`cdl_conformal_distance` in `null_rays.py` is $X(\xi) = \int d\xi/\rho$ through the wall, and `--verify` holds the wall chart's rays to $\psi \pm X$.

The conformal diagram uses the same $X$: $U = -e^{X - \psi}$ and $V = e^{X + \psi}$ are null outside the light cone and continue inside it as $e^{\eta \mp \chi}$, and $p = \arctan(U/s)$, $q = \arctan(V/s)$ draw both sides of the wall with one map.
For the decay into flat space $s = 8/5$ makes the outside de Sitter's own square, and for the decay of flat space $s = 1$ puts the wall at the same place, $X = 2\arctan(1/2)$, at the moment it appears.
The docstring of `coleman_de_luccia` in `conformal.py` has the details and its checks.
Only the half of each diagram after the bubble appears is drawn, since the earlier half of the classical solution, a bubble contracting to rest, is replaced by the tunnelling.

The embedding diagram is the equator of a moment of the open universe inside the anti-de Sitter bubble, the hyperbolic plane of radius $a = \ell\sin(c\tau/\ell)$, drawn in three dimensional Minkowski space as Milne's is, as a movie from $c\tau = \pi\ell/6$ to $5\pi\ell/6$.
Those moments are marked on the drawings of that bubble that hold them; `HIDDEN` in `slices.py` names the others.
No figure in three dimensions is drawn, since every cone of this spacetime lies in a plane of two coordinates already drawn.
