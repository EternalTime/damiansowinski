# Murphy's universe with viscosity

George L. Murphy, "Big-Bang Model Without Singularities", Phys. Rev. D 8, 4231 (1973).
The full text was not read (closed at APS, ADS and every open route); the abstract was.
The model is taken as Brevik and Grøn restate it (Astrophys. Space Sci. 2013, arXiv:1306.5634, section 5): a flat Friedmann universe of a fluid with $p = w\rho c^2$ and a bulk viscosity $\zeta = \beta\rho c^2$ proportional to its energy density, with $\beta$ a constant time.
Belinskii and Khalatnikov (Sov. Phys. JETP 42, 205 (1975)) name the same solution the isotropic separatrix of their flat Bianchi I system, the Friedmann parabola with second viscosity alone.

## The integral

With $\epsilon = \rho c^2$ and Eckart's viscous pressure $\Pi = -3\zeta H$, Friedmann's equation and the conservation of energy are

$$H^2 = \frac{8\pi G}{3c^2}\epsilon, \qquad \dot\epsilon = -3H\left(\epsilon + p - 3\zeta H\right) = -3H\epsilon\left(1 + w - 3\beta H\right).$$

The first makes $\dot\epsilon/\epsilon = 2\dot H/H$, so

$$\dot H = -\frac{3}{2}H^2\left(1 + w - 3\beta H\right),$$

which is Brevik and Grøn's (41).
It has the fixed point $H_s = (1 + w)/3\beta$, the steady state.
Write $q = 3(1 + w)/2$, so $H_s = 2q/9\beta$ and $\dot H = -qH^2(1 - H/H_s)$.

Along the solution $da/a = H\,dt = H\,dH/\dot H$, so

$$\frac{da}{a} = -\frac{dH}{qH\left(1 - H/H_s\right)}, \qquad \ln a = \frac{1}{q}\ln\frac{H_s - H}{H} + \text{const},$$

and with the constant chosen so that $a = 1$ where $H = H_s/2$,

$$a^q = \frac{H_s}{H} - 1, \qquad H = \frac{H_s}{1 + a^q}.$$

Brevik and Grøn's (44) prints the same relation with the opposite sign, $a^{3(1+w)/2} = 3\gamma - (1 + w)/H$ in their units, which makes $a$ fall as $H$ falls; the form here is the one the field equations hold (murphy_check in print_charts.py).

The proper time follows from $dt = da/(aH) = (1 + a^q)\,da/(H_sa)$:

$$H_st = \ln a + \frac{a^q}{q},$$

with $t = 0$ at the one moment where $\ln a = -a^q/q$.
Multiplying by $q$ and exponentiating, $a^qe^{a^q} = e^{qH_st}$, so $a^q = W(e^{qH_st})$ with $W$ Lambert's function, the cosmic chart's scale factor in closed form; the diagrams of that chart draw it so.
As $t \to -\infty$, $a \to 0$ like $e^{H_st}$ and $H \to H_s$: the steady state, de Sitter's expansion in its flat chart.
As $t \to \infty$, $a^q \to qH_st$, so $a \propto t^{2/3(1+w)}$ and $H \to 2/3(1 + w)t$, the Friedmann expansion of the same fluid without viscosity.

## The charts

The scale factor is monotonic, so it serves as the time:

$$ds^2 = -\frac{c^2\left(1 + a^q\right)^2}{H_s^2a^2}da^2 + a^2\left(dx^2 + dy^2 + dz^2\right),$$

in Cartesian and in spherical comoving coordinates.
The cosmic chart is the proper time with $a(t)$ left a function, as Murphy and Brevik and Grøn write it, held to $a' = (H_s/c)\,a/(1 + a^q)$ by the check.

murphy_check holds each chart, at six random points with random $q$, $H_s$ and $c$, to the viscous fluid at rest: a diagonal mixed Einstein tensor, $G^0{}_0 = -3H^2/c^2$, and $G^i{}_i = -\left(w - (1 + w)H/H_s\right)G^0{}_0$ with $w = 2q/3 - 1$.

## Curvature

$$R = \frac{6H_s^2\left(2 + (2 - q)a^q\right)}{c^2\left(1 + a^q\right)^3}, \qquad K = \frac{12H_s^4\left(\left(1 + a^q\right)^2 + \left(1 + (1 - q)a^q\right)^2\right)}{c^4\left(1 + a^q\right)^6}.$$

Both are finite for every $a > 0$ and tend to de Sitter's $12H_s^2/c^2$ and $24H_s^4/c^4$ as $a \to 0$; the Weyl tensor vanishes, as for every Friedmann universe.

## The past

A radial light ray has affine length $\lambda \propto \int a\,dt = \int (1 + a^q)\,da/H_s$, which is finite as $a \to 0$.
So light rays reach $t \to -\infty$ at a finite affine distance, as Borde, Guth and Vilenkin's bound requires of any universe whose Hubble rate stays positive (Phys. Rev. Lett. 90, 151301 (2003), their (3)-(5)): the past is null incomplete though the curvature is finite there.

The conformal time $c\,d\eta = c\,dt/a$ runs over $(-\infty, \infty)$ when $q > 1$, that is $w > -1/3$: $c\eta = (c/H_s)\left(a^{q-1}/(q - 1) - 1/a\right)$.
So the conformal diagram is the whole diamond of Minkowski space, its past null boundary at $a = 0$.

## The embedding diagram

Every moment is flat, so a moment embeds as a plane and shows nothing of the expansion.
The diagram stacks a ring of comoving observers at $r = c/H_s$ in the plane $y = 0$ up the axis of the proper time, from $H_st = -3$ to $3$ for radiation: the ring's proper radius $ac/H_s$ at the height of its time.
With $a^2 = W(e^{2H_st})$ the world tube narrows as $e^{H_st}$ toward the steady state and widens as $t^{1/2}$ late.
