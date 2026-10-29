# Ozsváth's homogeneous universes of dust

The four charts of `ozsvath_dust.json`, where each comes from and how each is checked.
`print_charts.py --metric ozsvath_dust` prints them and `verify_metrics.py --system ozsvath_dust/<chart>` checks them.

## Step 1. The sources and what was read

Ozsváth's papers of 1965 and 1970 and Farnsworth and Kerr's of 1966 are behind the publisher's wall and were not read; their records were checked on Crossref.
The line elements come from papers that print them and were read in full:
Sáez, Mengual and Ferrando, CQG 41, 205013 (2024), eq. (48), for classes II and III;
Wylleman, J. Phys. Conf. Ser. 229, 012073 (2010), eq. (3), for Farnsworth and Kerr's form of class III;
Costa and Natário, JMP 46, 082501 (2005), Theorem 4.3(3), for class I.
Sáez, Mengual and Ferrando show that the book's printing of classes II and III, their eq. (47), is no perfect fluid, and print the corrected eq. (48), which is used here.

Every chart is a metric left invariant on a four-dimensional group, so its curvature in the invariant frame is constant.
The checks below were first made in that frame from the structure constants alone, which takes a second, and then in the coordinates by the printer and the checker.

## Step 2. Dust and a cosmological constant

With $\kappa = 8\pi G/c^2$, dust of density $\rho$ and a cosmological constant $\Lambda$ give $R^\mu{}_\nu = \kappa\rho\left(u^\mu u_\nu + \tfrac12\delta^\mu_\nu\right) + \Lambda\delta^\mu_\nu$.
So $R^\mu{}_\nu$ has one eigenvalue $\Lambda - \kappa\rho/2$ on the timelike $u$ and a triple one $\Lambda + \kappa\rho/2$, and $R^\mu{}_\nu - (\Lambda + \kappa\rho/2)\delta^\mu_\nu$ has rank one and trace $-\kappa\rho$.
`ozsvath_dust_check` in `print_charts.py` holds each chart to that, with these values:

| chart | $\Lambda$ | $\kappa\rho$ | positive density |
| --- | --- | --- | --- |
| class I | $1/2L^2(1 + \cos^2\beta)$ | $\cos2\beta/L^2(1 + \cos^2\beta)$ | $\lvert\beta\rvert < \pi/4$ |
| class II | $-1/2a^2(1 + \cos^2\gamma)$ | $-\cos2\gamma/a^2(1 + \cos^2\gamma)$ | $\pi/4 < \lvert\gamma\rvert < \pi/2$ |
| class III, Farnsworth-Kerr | $-1/2a^2$ | $1/a^2$ | every $s$, every $v$ |

In Sáez, Mengual and Ferrando's $s = \sin\gamma/\sqrt2$, class II has $\kappa\rho\,a^2 = (4s^2 - 1)/2(1 - s^2)$ and $\Lambda a^2 = -1/4(1 - s^2)$, positive density exactly on their range $1/4 < s^2 < 1/2$.
At $s^2 = 1/2$ their $\lambda = \sqrt{1 - 2s^2}$ vanishes and the values are class III's, so the two classes meet there.
Class II's dust is tilted: its velocity lies in the plane of $dx - \sinh y\,dz$ and $du + \lambda(dx - \sinh y\,dz)$, at $\sqrt{14}/5$ of the speed of light along the second at $s = 3/5$.

## Step 3. Class I, built from Costa and Natário

Costa and Natário write a stationary metric $-(c\,dt - A)^2 + \gamma$ on a group, with $\gamma$ left invariant and an orthonormal frame $X_i$ whose brackets are $[X_2, X_3] = aX_1$ and its cyclic companions with $b$ and $c$.
Their class I has $c = a - b$, $a/4 < b < 3a/4$, a vanishing gravitational field and a gravitomagnetic field $H$ along $X_1$ with $H^2 = 4b(a - b) - a^2/2$, a perfect fluid with $p = -a^2/64\pi$ and $\rho = (32b(a - b) - 5a^2)/64\pi$.
With the Euler forms of $SU(2)$, $d\sigma_1 = -\sigma_2\wedge\sigma_3$ and its cyclic companions, the frame is $X_i^\flat = k_i\sigma_i$ with $k_1 = 1/\sqrt{bc}$, $k_2 = 1/\sqrt{ac}$, $k_3 = 1/\sqrt{ab}$, and $A$ is $H/a$ times $k_1\sigma_1$.
The chart puts the axis on $\sigma_3 = d\psi + \cos\theta\,d\phi$, sets $b/a = (1 + \epsilon)/2$, and writes $\epsilon = \sin\beta/\sqrt2$ and $L^2 = 2/a^2(1 - \epsilon^2)$, so that the twist $Lh$, with $h^2 = (1 - 2\epsilon^2)/(1 - \epsilon^2)$ in units of $R = L\sqrt{1 - \epsilon^2}$, becomes $L\cos\beta$.
That relation between the twist and the anisotropy is an elliptic curve, so no rational parameter removes the root, and a root whose radicand leads with a negative term is read on its other branch; the angle $\beta$ removes it.
The fluid of Costa and Natário reads as dust with $\Lambda = -8\pi Gp/c^2$, which is the table's: $\Lambda = a^2/8$ and $\kappa\rho = a^2(4b(a - b)/a^2 - 3/4)$ in their units.

## Step 4. The members of more symmetry

Class III at $s = 0$ is $a^2\left(-2(dx - \sinh y\,dz)^2 + dy^2 + \cosh^2y\,dz^2 + du^2\right)$, Gödel's universe in the form whose space of $(x, y, z)$ is the group $SL(2, \mathbb{R})$ stretched along its timelike direction, with $\kappa\rho = 1/a^2$ and $\Lambda = -\kappa\rho/2$.
Krasiński (JMP 39, 380, 1998, sec. 9) finds Gödel's universe as the shear-free member of Ozsváth's class III, and Wylleman as $v = 0$ of eq. (3).
Class I at $\beta = 0$ has a vanishing Weyl tensor, Ricci scalar $3/2L^2$ and Kretschmann scalar $3/4L^4$, those of Einstein's static universe of radius $2L$, and $\Lambda = 1/4L^2 = \kappa\rho/2$, Einstein's relation.
It is Costa and Natário's point $b = a/2$, where their Einstein family $c = b$ and their class I $c = a - b$ meet.

## Step 5. Farnsworth and Kerr's form is class III

Wylleman's eq. (3) is $A^{-2}$ times $du^2 - 2(dt + e^x dy)^2 + (1 + \tanh v)(\cos t\,dx + \sin t\,e^x dy)^2 + (1 - \tanh v)(-\sin t\,dx + \cos t\,e^x dy)^2$, since $e^{\pm v}/\cosh v = 1 \pm \tanh v$.
Its four forms have the structure constants of class III's four with every sign reversed, so negating each frame vector carries one to the other, and the frame metrics $\mathrm{diag}(-2, 1 + s, 1 - s, 1)$ agree at $s = \tanh v$ and $a = 1/A$.
Two metrics left invariant on groups whose frames have the same structure constants and the same frame metric are locally isometric.
Both charts print $R = -1/a^2$ and $K = 3/a^4$ at every $s$.

## Step 6. Printing

Class III prints in $\sinh y$ and $\cosh y$ through `chart_printer.hyperbolic`, and Farnsworth-Kerr in $\sinh v$ and $\cosh v$.
Class II is printed in Sáez, Mengual and Ferrando's $s$ and $\lambda$: `ozsvath_class_two_printing` writes $\sin\gamma = \sqrt2 s$ and $\cos\gamma = \lambda$ and reduces by $\lambda^2 = 1 - 2s^2$; the reader reads $s$ and $\lambda$ back by their definitions, so every printed value is still compared with the one it came from.
