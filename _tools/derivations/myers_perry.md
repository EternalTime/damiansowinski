# The Myers-Perry black hole

Myers and Perry's rotating black hole has one spin for each independent plane of rotation, $\lfloor (D-1)/2 \rfloor$ of them in $D$ spacetime dimensions.
A chart needs its coordinates named one by one, so five charts are published, four in $D = 5$ and one in $D = 6$, and the prose carries the general $D$.

## Step 1: the parameters and the sign of the spin

The mass parameter is $\mu = 16\pi GM/((D-2)\Omega_{D-2}c^2)$ and each spin parameter is $a_i = (D-2)J_i/(2Mc)$, Emparan and Reall's (34) with $c$ restored.
$\mu$ carries the dimension $L^{D-3}$: $8GM/(3\pi c^2)$, an area, in five dimensions and $3GM/(2\pi c^2)$, a volume, in six, and `DIMENSIONS` declares it so chart by chart.
Myers and Perry's own metric has $(dt + a\sin^2\theta\,d\phi)^2$, rotation toward decreasing $\phi$.
Emparan and Reall write $(dt - a\sin^2\theta\,d\phi)^2$, rotation toward increasing $\phi$, which is Kerr's sign in its own chart here, and every chart takes it.

## Step 2: one spin, five dimensions and six

`boyer_lindquist` and `boyer_lindquist_six` are Emparan and Reall's (32) at $d = 5$ and $d = 6$, with $\Sigma = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 + a^2 - \mu/r^{D-5}$.
$\theta$ runs from $0$ to $\pi/2$: $\theta = \pi/2$ is the plane of rotation, where the transverse sphere has shrunk to a point, and $\theta = 0$ is the space transverse to the rotation, where the circle of $\phi$ has.
In five dimensions the horizon is $r_+ = \sqrt{\mu - a^2}$; in six it is the one positive root of $r^3 + a^2r = \mu$, which exists for every $a$.
The six-dimensional chart is checked to be Tangherlini's on its plane of $t$ and $r$ at $a = 0$ and flat at $\mu = 0$.

## Step 3: the ingoing chart

Myers's Eddington-like coordinates, his (1.27) for odd $d$, are $dt = dt_+ - \mu r^2\,dr/(\Pi - \mu r^2)$ and $d\phi = d\phi_+ + \Pi\,a\,dr/((\Pi - \mu r^2)(r^2 + a^2))$ with $\Pi = r^2(r^2 + a^2)$ for one spin.
`ingoing_kerr` takes the null $v = ct_+ + r$ for its time, so $dv = c\,dt + (r^2 + a^2)\,dr/\Delta$ and $d\tilde\phi = d\phi + a\,dr/\Delta$, as the ingoing Kerr coordinates do in four dimensions.
The metric is the flat one in spheroidal coordinates plus $\mu/\Sigma$ times the square of $dv - a\sin^2\theta\,d\tilde\phi$, the Kerr-Schild form, and `myers_perry_check` pulls the Boyer-Lindquist chart back onto it slot by slot.

## Step 4: two spins and equal spins

`two_spins` is Myers's (1.66) with $\phi_i \to -\phi_i$: $\Sigma = r^2 + a^2\cos^2\theta + b^2\sin^2\theta$ and $g_{rr} = r^2\Sigma/((r^2 + a^2)(r^2 + b^2) - \mu r^2)$, checked to be the one-spin chart at $b = 0$.
At $b = a$ the potential $\sin^2\theta\,d\phi + \cos^2\theta\,d\psi$ is $\tfrac12(d\Psi + \cos\Theta\,d\Phi)$ with $\Theta = 2\theta$, $\Phi = \psi - \phi$ and $\Psi = \psi + \phi$, the Euler angles of the 3-sphere, and $d\theta^2 + \sin^2\theta\,d\phi^2 + \cos^2\theta\,d\psi^2$ is a quarter of its round metric.
With $\rho^2 = r^2 + a^2$, $r^2dr^2 = \rho^2d\rho^2$ and $g_{\rho\rho} = (1 - \mu/\rho^2 + \mu a^2/\rho^4)^{-1}$: `equal_spins`, which is Murata and Soda's (1) with their $2\mu$ written $\mu$ and the sign of $a$ reversed.
The chart reaches below $r = 0$ of the Boyer-Lindquist radius, to the singularity $\rho = 0$, and its horizons are $2\rho_\pm^2 = \mu \pm \sqrt{\mu^2 - 4\mu a^2}$.
`myers_perry_check` pulls the two-spin chart at $b = a$ back onto it.

## Step 5: the curvature and the printing

Every chart is Ricci flat, which the script checks before it writes, so the Weyl tensor is the Riemann tensor.
The Kretschmann scalar is Myers's (1.69), $24\mu^2(4r^2 - 3\Sigma)(4r^2 - \Sigma)/\Sigma^6$ in five dimensions, printed with its factors written out, and $24\mu^2(\rho^2 - 4a^2)(3\rho^2 - 4a^2)/\rho^{12}$ for equal spins.
`myers_perry_pretty` factors each value with $\sin^2\theta$ written as $1 - \cos^2\theta$, the angle $\Sigma$ is written in, so that $\Sigma$ and the two-spin $\Delta$ come out as factors and keep the order the line element writes them in.

## Step 6: the diagrams

On the plane of $t$ and $r$ at $\theta = 0$ with one spin the metric is $-f\,c^2dt^2 + dr^2/f$ with $f = \Delta/(r^2 + a^2)$, and no Christoffel symbol turns a ray out of it.
In five dimensions $1/f = 1 + \mu/(r^2 - r_+^2)$, so $r_* = r + (\mu/2r_+)\ln|(r - r_+)/(r + r_+)|$, which vanishes at $r = 0$, and the surface gravity is $r_+/\mu$: Schwarzschild's square, as Myers finds for one vanishing spin in odd dimensions with $\mu > a^2$.
The surface $r = 0$ is the locus of a conical singularity with the Kretschmann scalar $72\mu^2/a^8$ on it, and the curvature diverges on its rim $\theta = \pi/2$.
In six, $1/f = 1 + \mu/(r^3 + a^2r - \mu)$ and $r_* = r + \mathrm{Re}\sum_i A_i\ln(1 - r/r_i)$ over the three roots of the cubic, $A_i = \mu/(3r_i^2 + a^2)$.

On the plane of rotation the circles of $\phi$ are divided out, as for the BTZ hole: $dr/d(ct) = \pm\Delta/\sqrt{r^4 + a^2r^2 + \mu a^2}$ in five dimensions.
For equal spins the Hopf fibre $\psi$ is divided out, and $d\rho_*/d\rho = \rho^2\sqrt{\rho^4 + \mu a^2}/(\rho^4 - \mu\rho^2 + \mu a^2)$, an elliptic integral, which `MyersEqual` takes as Kerr's `EquatorTower` does: the residues at the four roots $\pm\rho_\pm$, which sum to zero, and a smooth remainder integrated numerically.
The two-spin chart has no spacetime diagram: on either plane of rotation the other plane's circle has shrunk to a point, where its angle is no coordinate and the lift of a ray of no angular momentum misses the geodesic equation, and off both planes a ray moves in $\theta$.
The six-dimensional plane of rotation is drawn at $a = 1.5\,\mu^{1/3}$ on a box of $4\,\mu^{1/3}$: at $a = 2\,\mu^{1/3}$ on a box of $3$ the generator's central differences, of step $10^{-5}$, miss the geodesic equation by $1.7 \times 10^{-4}$ at its sample nearest $r = 0$, where the components vary as $r^{-3}$.
Figueras, Kunesch, Lehner and Tunyasuvunakool evolved $1.5 \le a/\mu^{1/3} \le 2$, and the first unstable mode sets in at $1.572$.

The embedding diagram has four views.
In five dimensions with one spin, the plane of rotation at constant $t$ has the circumference radius $\sqrt{r^2 + a^2 + \mu a^2/r^2}$, $\mu/r_+$ at the throat, and the transverse plane is $z = \sqrt{\mu}\,\mathrm{arcosh}(r/r_+)$ with a throat of radius $r_+$.
For equal spins the surface of $\rho$ and the fibre at $\theta = \pi$ sweeps $\psi = 3\phi$, since $d\psi - d\phi$ must go once round a fibre of period $4\pi$, and the fibre's radius is $\sqrt{\rho^4 + \mu a^2}/\rho$, $\sqrt{\mu}$ at the throat.
In six dimensions the transverse plane has $dz/dr = \sqrt{\mu/(r^3 + a^2r - \mu)}$ and a finite height.
The plane of rotation in six dimensions has no embedding in one piece at this spin: between the throat and $r \approx \mu^{1/3}$ its circles shrink outward faster than the distance out to them, as at $r = 0.6\,\mu^{1/3}$, where $(d\rho/dr)^2 = 50$ and $g_{rr} = 0.38$, and far out they grow more slowly, so neither flat space nor Minkowski space carries the whole slice.
