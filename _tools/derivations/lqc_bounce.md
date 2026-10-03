# The bounce of loop quantum cosmology

The flat Friedmann universe of a massless scalar field in the effective dynamics of loop quantum cosmology, in three charts.
`print_charts.py` writes every chart's mathematics and `lqc_bounce_check` holds each to the equations below.

## Step 1: the modified Friedmann equation

Ashtekar, Pawlowski and Singh, Phys. Rev. D 74, 084003 (2006), Appendix B, derive from the effective Hamiltonian of the improved dynamics
$$H^2 = \frac{8\pi G}{3}\rho\left(1 - \frac{\rho}{\rho_c}\right),$$
with $\rho = p_\phi^2/2|p|^3$ the density of the scalar field.
They print $\rho_c = \sqrt{3}/16\pi^2\gamma^3G^2\hbar \approx 0.82\rho_{\rm Pl}$.
Ashtekar and Wilson-Ewing, Phys. Rev. D 79, 083535 (2009), take for a homogeneous universe the area gap $4\sqrt{3}\pi\gamma\ell_p^2$, twice the smallest eigenvalue of the area, which halves the critical density to $\sqrt{3}/32\pi^2\gamma^3G^2\hbar \approx 0.41\rho_{\rm Pl}$, the value of Ashtekar and Singh's status report, Class. Quantum Grav. 28, 213001 (2011), and of the 2010 arXiv revision of Ashtekar, Corichi and Singh.
The page keeps $\rho_c$ symbolic and quotes 0.41.

A massless scalar field has $\rho \propto a^{-6}$, so with $a = 1$ at the bounce $\rho = \rho_c/a^6$.
Write $t_b = 1/\sqrt{24\pi G\rho_c}$; then $8\pi G\rho_c/3 = 1/9t_b^2$ and
$$\left(\frac{\dot a}{a}\right)^2 = \frac{1}{9t_b^2a^6}\left(1 - \frac{1}{a^6}\right).$$

## Step 2: the cosmic chart

With $u = a^6$ the equation is $\dot u^2 = 4(u - 1)/t_b^2$, so $u = 1 + t^2/t_b^2$ with the bounce at $t = 0$:
$$a(t) = \left(1 + \frac{t^2}{t_b^2}\right)^{1/6} = \left(1 + 24\pi G\rho_c t^2\right)^{1/6}.$$
Cai and Wilson-Ewing, JCAP 03 (2014) 026, print it in units $c = \hbar = 1$ with $M_{\rm Pl}^2 = 1/8\pi G$: after the bounce, $(3\rho_ct^2/M_{\rm Pl}^2 + 1)^{1/6}$, and their ekpyrotic law at $p = 1/3$, the same, before it.
`lqc_bounce_check` confirms that it solves the equation of Step 1 for every $t$.
The chart is $ds^2 = -c^2dt^2 + a^2(dx^2 + dy^2 + dz^2)$, the form of Ashtekar and Singh's harmonic metric with unit lapse, and the comoving spherical chart is the same with $dr^2 + r^2d\Omega^2$.

Every component is $a^{2j}$ times a rational function of $t$: `lqc_scale_factor` picks the power that leaves $t_b^2 + t^2$ whole.
The Ricci scalar is $R = 2(3t_b^2 - t^2)/3c^2(t_b^2 + t^2)^2$ and the Kretschmann scalar $4(9t_b^4 - 12t_b^2t^2 + 5t^4)/27c^4(t_b^2 + t^2)^4$, both finite, the latter greatest at the bounce, $4/3c^4t_b^4$.
The Weyl tensor vanishes.

## Step 3: the harmonic chart

Ashtekar, Corichi and Singh, Phys. Rev. D 77, 024046 (2008), take the scalar field as time already in the classical theory and find the volume $V_{\rm min}\cosh(\sqrt{12\pi G}\,\phi)$ about the bounce for every state, and the improved dynamics' effective equation $dv/d\phi = \sqrt{12\pi G}(1 - \rho/\rho_c)^{1/2}v$ integrates to the same.
Ashtekar and Singh's status report writes the metric in the harmonic time $\tau$, $\Box\tau = 0$, as $-a^6d\tau^2 + a^2d\vec x^2$, lapse $a^3$.
With $c\,dt = a^3c\,d\tau$ and $a^6 = 1 + t^2/t_b^2$, $t = t_b\sinh(\tau/t_b)$ and $a = \cosh^{1/3}(\tau/t_b)$, and the volume $a^3 = \cosh(\tau/t_b)$ is Ashtekar, Corichi and Singh's law with $\tau/t_b = \sqrt{12\pi G}(\phi - \phi_b)$.
`lqc_bounce_check` holds the chart to $\Box\tau = 0$ and to the cosmic chart pulled back; `slices.py` checks the same pullback numerically.

## Step 4: the drawings

Light obeys $dx = \pm c\,dt/a$, so $x \pm \eta$ is constant along a ray with $\eta = \int_0^t dt/a = t\,{}_2F_1(1/6, 1/2; 3/2; -t^2/t_b^2)$, which grows as $(3/2)t_b^{1/3}|t|^{2/3}$ without limit both ways.
`null_rays.py --verify` holds every ray of the three spacetime diagrams to it, with $t = \sinh\tau$ in the harmonic chart.
Since $\eta$ covers the whole real line, $p, q = \arctan((\eta \mp r)/ct_b)$ map the comoving chart onto the whole of Minkowski's half diamond, which `conformal.py` draws.
Each moment of constant $t$ is flat, so the embedding diagram is the equatorial plane of the comoving chart at five moments and a movie through the bounce, with the comoving observer at $r$ drawn at radius $ar$.
The Hubble sphere of the comoving diagram is $r = 3c(t_b^2 + t^2)/a|t|$, nearest the centre at $t = \pm\sqrt{3/2}\,t_b$, $r \approx 5.26\,ct_b$.
