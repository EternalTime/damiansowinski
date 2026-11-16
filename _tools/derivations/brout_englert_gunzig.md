# The causal universe of Brout, Englert and Gunzig

## What is stated, and from which source

The metric is conformal to the flat one, $g_{\mu\nu} = (\kappa\phi^2/6)\,\eta_{\mu\nu}$, with a factor that depends on $\tau^2 = t^2 - r^2/c^2$ alone inside the forward light cone of one event and is constant outside it.
That is the essay of Brout, Englert and Gunzig, General Relativity and Gravitation 10, 1 (1979), doi:10.1007/BF00757018, read in the scan of the Gravity Research Foundation.
We write the square root of the factor as $\Omega(\tau)$, so that $ds^2 = \Omega^2\left(-c^2dt^2 + dr^2 + r^2d\theta^2 + r^2\sin^2\theta\,d\varphi^2\right)$, with $\Omega = 1$ outside the cone.
The azimuth is $\varphi$ throughout, since $\phi$ is the authors' scalar field.

The fireball is taken with the correction $\delta$ to the density of created particles set to zero.
The essay states that the solution is then "an open universe which closely approximates a de Sitter space", with the factor unbounded at a finite $\tau_\infty$.
Gott, IAU Symposium 104, 475 (1983), doi:10.1017/S0074180900039565, states the scale factor $a \propto \sinh(t/t_{ex})$ and the open geometry $k = -1$ inside the cone.
Spindel and Vandepeutte, Liège International Astrophysical Colloquium 26, 121 (1986), state $H = 2/\tau_\infty$ on their page 128 and write the open line element in conformal time as their equation (18).
The article in Annals of Physics 115, 78 (1978) was not read, and nothing here rests on its formulas.

The free expansion is the radiation era alone: an open Friedmann universe with pressure one third of its energy density, which Gott names as the phase that follows the fireball.
The essay states that the two stages join smoothly and that continuity fixes the joining.
The universe of 1978 also holds massive matter after the fireball; the charts here stop at the radiation era.

## The three coordinate systems

With $t = \tau\cosh\chi$ and $r = c\tau\sinh\chi$ the flat metric is $-c^2d\tau^2 + c^2\tau^2\left(d\chi^2 + \sinh^2\chi\,d\Omega^2\right)$, so the scale factor of the hyperboloids is $a = c\tau\Omega$, the proper time of the comoving world lines obeys $dT = \Omega\,d\tau$, and the conformal time is $\eta = \ln(\tau/\tau_\infty)$.

The fireball: $\Omega = \left(1 - \tau^2/\tau_\infty^2\right)^{-1}$.
Then $T = \tau_\infty\,\mathrm{artanh}(\tau/\tau_\infty)$, $a = (c/H)\sinh HT$ with $H = 2/\tau_\infty$, and $a = -c/(H\sinh\eta)$ on $\eta < 0$.
Its Ricci tensor is $(3H^2/c^2)\,g_{\mu\nu}$, de Sitter space in hyperbolic slices.

The free expansion: $\Omega = \Omega_f\left(1 - \tau_b^2/\tau^2\right)$.
Then $a = c\sqrt{(T - T_b)^2 - b^2}$ and $a = bc\sinh(\eta - \eta_b)$, with $b = 2\tau_b\Omega_f$ and $\eta_b = \ln(\tau_b/\tau_\infty)$.
Its Ricci scalar vanishes and its Einstein tensor is that of a fluid at rest on the hyperboloids with $G_{TT} = 3b^2c^2/a^4$ in the chart $x^0 = cT$.

## The joining

Let $x = \tau_1/\tau_\infty$ and $\Omega_1 = 1/(1 - x^2)$.
The fireball has $d\Omega/d\tau = 2\tau\Omega^2/\tau_\infty^2$ and the free expansion $d\Omega/d\tau = 2\Omega_f\tau_b^2/\tau^3$.
Equating $\Omega$ and $d\Omega/d\tau$ at $\tau_1$ gives $\Omega_f\tau_b^2 = x^2\Omega_1^2\tau_1^2$ and $\Omega_f = \Omega_1 + x^2\Omega_1^2 = \Omega_1^2$, so

$$\Omega_f = \left(1 - \tau_1^2/\tau_\infty^2\right)^{-2}, \qquad \tau_b = \tau_1^2/\tau_\infty.$$

In conformal time that is $\eta_b = 2\eta_1$ and $b = 1/(H\sinh^2\eta_1)$.
In proper time, with $\tanh(HT_1/2) = x$, it is $b = \sinh^2(HT_1)/H$ and $T_b = T_1 - \sinh(2HT_1)/(2H)$.
The authors' figure marks $e^\lambda = \Omega^2 = 10^{16}$ where creation ends and $10^{32}$ at late times, which is $\Omega_f = \Omega_1^2$.

## The light cone

On the cone $\Omega = 1$ on both sides.
Inside it $\Omega = 1 + \tau^2/\tau_\infty^2 + \dots$, so the derivative of $\Omega$ along $\tau^2$ is $1/\tau_\infty^2$ inside and zero outside.
The metric is continuous across the cone and its first derivative is not, in the form with $\delta = 0$.
The essay does not treat the interval $\tau \lesssim \hbar/mc^2$.

## The drawings

Every drawing takes $\tau_1 = 0.8\,\tau_\infty$, which keeps the two stages on one figure; the authors' own value for $m = 1$ GeV has $\Omega_1 = 10^8$.
Each is drawn in units of $1/H$ and $c/H$, in which $\tau_\infty = 2$, $\tau_1 = 1.6$, and $HT_1 = 2\ln 3$.

Null rays are straight lines at 45° on the flat charts and on the conformal time charts, and on the proper time charts they conserve $\eta \pm \chi$ with $\eta = \ln\tanh(HT/2)$ in the fireball and $\cosh(\eta - \eta_b) = (T - T_b)/b$ after it.

The conformal diagram is flat space's triangle, $p, q = \arctan(H(t \mp r/c))$.
The factor of the free expansion tends to the constant $\Omega_f$, so its future infinity is flat space's.
Every hyperboloid of constant $\tau$ runs from the central world line to the point where the light cone of the origin meets future null infinity, as in the diagram of Gott's Figure 1.

The embedded surface is the plane $\theta = \pi/2$ at one time $t$ of the flat background, with circles of radius $\rho = \Omega r$.
A surface of revolution in flat space carries it where $g_{rr} \ge (d\rho/dr)^2$.
In the fireball $g_{rr} - (d\rho/dr)^2 = 4r^2\Omega^4\left(\tau_\infty^2 - t^2\right)/\tau_\infty^4$, and in the free expansion it is $4r^2\Omega_f^2\tau_b^2\left(\tau^4 - t^2\tau_b^2\right)/\tau^8$, which is positive on $\tau \ge \tau_1$ exactly when $t \le \tau_\infty$.
So the plane is embedded for $t \le \tau_\infty$, and the drawing takes $t = 0.9\,\tau_\infty$, later than the end of creation on the central world line.
At the cone $d\rho/dr = 1 - 2t^2/\tau_\infty^2$ inside and $1$ outside, each over $\sqrt{g_{rr}} = 1$, so the surface has a crease there.
The hyperboloids of constant $\tau$ have $g_{\chi\chi} < (d\rho/d\chi)^2$ at every $\chi > 0$, so no surface of revolution in flat space carries them.

The curvature figure plots the Kretschmann scalar and the Ricci scalar along the proper time.
In the fireball $K = 24H^4/c^4$ and $R = 12H^2/c^2$; in the free expansion $R = 0$ and $K = 24b^4c^4/a^8$, which is $24H^4/c^4$ at $T_1$, since $a_1 = (c/H)\sinh HT_1$ and $b = \sinh^2(HT_1)/H$.
