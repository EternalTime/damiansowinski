# Mukhanov and Brandenberger's universe with limiting curvature

The theory is Viatcheslav Mukhanov and Robert Brandenberger's of 1992, Phys. Rev. Lett. 68, 1969, and its universes are worked out by Robert Brandenberger, Viatcheslav Mukhanov and Andrew Sornborger in Phys. Rev. D 48, 1629 (1993), arXiv:gr-qc/9303001, whose equation numbers are used below.
No universe of the theory is known in closed form, so the page carries the closed Robertson-Walker chart with its scale factor left free and the closed slicing of de Sitter space, the geometry every solution approaches, and draws a numerical solution that `limiting_curvature.py` beside this file integrates.

## The action and the sign of the root

In the site's signature $(-,+,+,+)$, with the standard contraction $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$, the action is $\sqrt{-g}\,[R + \phi_1I_1 - \phi_2\sqrt{I_2} - V_1(\phi_1) - V_2(\phi_2)]$, with $I_1 = R - \sqrt{3}\sqrt{I_2}$ and $I_2 = 4R_{\mu\nu}R^{\mu\nu} - R^2$.
On a Robertson-Walker universe of spatial curvature $k$, in units with $c = 1$, $I_2 = 12(\dot H - k/a^2)^2$, and the root that makes the theory work is the signed one, $\sqrt{I_2} = +\sqrt{12}(\dot H - k/a^2)$, which carries the sign of $\dot H - k/a^2$.
With it $I_1 = R - 6(\dot H - k/a^2) = 12(H^2 + k/a^2)$ identically on both branches, and varying $\phi_1$ and $\phi_2$ in the minisuperspace action gives the paper's constraints (2.25), $H^2 + k/a^2 = V_1'/12$, and (2.26), $\dot H - k/a^2 = -V_2'/\sqrt{12}$.
With the positive root neither comes out: at $k = 0$ it gives $I_1 = 12H^2 + 6\dot H - 6|\dot H|$, which is $12H^2$ only where $\dot H \ge 0$.
In the paper's signature $(+,-,-,-)$ the same root is $-\sqrt{12}(\dot H - k/a^2)$, since $R$ changes sign with the signature and $R_{\mu\nu}R^{\mu\nu}$ does not.
The paper never states the sign; its equations fix it.
`_tools/test_limiting_curvature.py` derives both constraints from the action with the signed root and shows the positive root gives neither, and the check of the comoving chart in `print_charts.py` holds the chart's own $I_1$ and $I_2$ to these forms at points where $\dot H - k/a^2$ takes each sign.

## The Hamiltonian equation

Varying the lapse $N$ and setting $N = 1$ gives
$$-\tfrac12(V_1 + V_2) + 3H^2(1 - 2\phi_1) + 3\frac{k}{a^2}(2\phi_1 + 1) = \sqrt{3}H\left(\dot\phi_2 + 3H\phi_2 - \frac{k\phi_2}{Ha^2}\right).$$
The paper prints the curvature term as $3(k/a^2)(4\phi_1 + 1)$ in (2.27) and as $3(k/a^2)(4\phi_1 - 1)$ in (4.2); each misprints one coefficient, and at $k = 0$ all three agree.
The variation is checked in `_tools/test_limiting_curvature.py`, and the solver uses the derived form.
The page's prose gives the two constraints, which the paper prints correctly, and none of the three forms of this equation.

## Integrating through the bounce

The Hamiltonian equation divides by $H$, which vanishes at a bounce.
At $H = 0$ the coefficient of $\dot\phi_2$ in it vanishes too, so it fixes $\phi_2$ there and leaves $\dot\phi_2$ free: with $\dot H_b = k/a^2 - V_2'/\sqrt{12}$ at the bounce, the equation linearised about the root reads $\dot\psi = \lambda\psi/(t - t_b)$ for the distance $\psi$ of $\phi_2$ from it, and $\lambda = (\sqrt{3}k/a^2 - V_2'/2)/(\sqrt{3}\dot H_b) = 1$ exactly, since $\sqrt{3}/\sqrt{12} = 1/2$.
So every solution passes through the bounce with $\phi_2$ finite, and the bounces of one radius form a family of one more number, the rate of $\phi_2$ there.
The solver steps over the division altogether: it integrates the variation of $a$, which is second order in $\phi_2$, has the coefficient $-2\sqrt{3}a^2$ in front of $\ddot\phi_2$ and divides by nothing, and keeps (2.25) and the Hamiltonian equation as checks, to $10^{-10}$ out to $t = 60/H_0$.
`limiting_curvature.py` writes that equation out, and the test derives it from the action and compares it at points.
Of the paper's (4.1), $\dot\phi_1 = -4\sqrt{3}HV_2'/V_1''$, the solver uses the form as printed.

## The potentials and the bounds

The potentials are the paper's (4.7) and (4.8), $V_1 = 12H_0^2\,\phi_1^2/(1 + \phi_1)\,(1 - \ln(1 + \phi_1)/(1 + \phi_1))$ and $V_2 = 2\sqrt{3}H_0^2\,\phi_2^2/(1 + \phi_2^2)$.
The slope of the first is $V_1'/12H_0^2 = \phi_1(\phi_1^2 + 2\phi_1 + 2 - 2\ln(1 + \phi_1))/(1 + \phi_1)^3$, which rises monotonically from 0 at $\phi_1 = 0$ toward 1 and never reaches it, so (2.25) keeps $H^2 + k/a^2$ below $H_0^2$, and a closed universe, $k = 1$, below which $a$ never falls, $a > 1/H_0$.
The largest size of the second's slope is $9H_0^2/4$, at $\phi_2 = \pm1/\sqrt{3}$, so (2.26) keeps $|\dot H - k/a^2|$ at most $9H_0^2/(4\sqrt{12}) = 3\sqrt{3}H_0^2/8$.
As $|\phi_2| \to \infty$ the slope $V_2'$ falls as $4\sqrt{3}/\phi_2^3$, and $\dot H - k/a^2 \to 0$, de Sitter space.

## The universe drawn

The drawings are of the closed universe, $k = 1$, without matter, through a bounce at $t = 0$ of radius $a_b = 2c/H_0$, in units $c = H_0 = 1$.
At the bounce $H = 0$, (2.25) gives $\phi_1 = 0.214787$, the Hamiltonian equation $\phi_2 = -0.705320$, the root nearest zero, and the rate of $\phi_2$ is set to zero, which makes the universe symmetric under $t \to -t$.
That symmetry takes $a$, $\phi_1$ and $\phi_2$ into themselves and $H$ and $\dot\phi_2$ into their negatives, so the solver integrates forward from the bounce and reads the past through it.
On either side $\phi_2$ runs off toward $-\infty$ and $\phi_1$ grows, as the paper describes the closed universe of its Fig. 10 and Fig. 11 near its de Sitter bounce, and the universe approaches de Sitter space: $\dot H - k/a^2$ is $0.030H_0^2$ at $t = 2/H_0$, $0.0027H_0^2$ at $10/H_0$ and $0.00047H_0^2$ at $60/H_0$, while $H$ climbs from $0.82H_0$ to $0.92H_0$.
$\phi_2$ stays negative, rising to $-0.359$ at $t = 0.538/H_0$ before it falls away, so $\phi_1$ grows monotonically, and with it $H^2 + c^2/a^2 = V_1'/12$, from $0.25H_0^2$ at the bounce to $0.850H_0^2$ at $t = 60/H_0$, below $H_0^2$ throughout.
The conformal time $\eta = \int c\,dt/a$ from the bounce tends to $0.8430$ as $t \to \infty$, so the whole history spans $1.686$, less than $\pi$, and no light ray crosses from a pole of the 3-sphere to its antipode.
