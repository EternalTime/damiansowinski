# Which radical carries the factor of three in the interior Schwarzschild solution

The interior Schwarzschild inverse metric had one wrong component, and we derive its correction.
Only that one component is at issue, not the whole spacetime.
That component did not invert the metric, and of the interior Schwarzschild disagreements listed in `_tools/derivations/audit-2026-09-18.md` it is the one that is real rather than a limit of the simplifier.

The line element carries the time term

$$-\frac{1}{4}\left(\sqrt{1 - \frac{r^2 r_s}{R^3}} - 3\sqrt{1 - \frac{r_s}{R}}\right)^{\!2}dt^2,$$

and the inverse metric carried

$$g^{tt} = -\frac{4}{\left(3\sqrt{1 - \frac{r_s}{R^3}r^2} - \sqrt{1 - \frac{r_s}{R}}\right)^{\!2}}.$$

The factor of three sits on the constant radical in the first and on the radius dependent radical in the second, so one of the two is wrong and the pair does not satisfy $g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu{}_\nu$.

We do not assume the line element wins.
Both readings are solutions of the same second order equation, and by Step 3 both of them join smoothly to the exterior Schwarzschild solution at $r = R$, which is why the error survived.
What separates them is the boundary condition at the surface of the star, imposed in Step 6, and three independent confirmations in Steps 8, 9 and 10 agree with the answer it gives.

The conclusion is that the line element is right and the inverse was wrong, and the corrected component is

$$g^{tt} = -\frac{4}{\left(3\sqrt{1 - \frac{r_s}{R}} - \sqrt{1 - \frac{r^2 r_s}{R^3}}\right)^{\!2}}.$$

The other fifty disagreements the checker reports for the interior Schwarzschild star are the same pair of radicals written as $\sqrt{(R-r_s)(R^3-r^2r_s)}$ where sympy writes $R^2\sqrt{1-r_s/R}\sqrt{1-r^2r_s/R^3}$, which sympy will not merge without being told the factors are positive.
They are not errors, and the correction leaves every one of them alone.

---

## Step 1. Conventions

The signature is $(-,+,+,+)$ and the chart coordinate is $x^0 = ct$, so every chart component of the metric is dimensionless and every factor of $c$ that a time index would carry is written out.
The connection is the Levi-Civita one,

$$\Gamma^\mu{}_{\nu\rho} = \tfrac{1}{2}g^{\mu\alpha}\left(\partial_\nu g_{\alpha\rho} + \partial_\rho g_{\alpha\nu} - \partial_\alpha g_{\nu\rho}\right),$$

the Riemann tensor is

$$R^\mu{}_{\nu\rho\sigma} = \partial_\rho \Gamma^\mu{}_{\nu\sigma} - \partial_\sigma \Gamma^\mu{}_{\nu\rho} + \Gamma^\mu{}_{\rho\lambda}\Gamma^\lambda{}_{\nu\sigma} - \Gamma^\mu{}_{\sigma\lambda}\Gamma^\lambda{}_{\nu\rho},$$

and the Ricci tensor is the standard contraction on the first lower index,

$$R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu},$$

which was settled for every spacetime on 18 September 2026.
The field equation is $G_{\mu\nu} = \dfrac{8\pi G}{c^4}T_{\mu\nu}$.

The source is a perfect fluid at rest in this chart,

$$T^\mu{}_\nu = \operatorname{diag}\left(-\rho c^2,\ p,\ p,\ p\right),$$

with $\rho$ constant, which is what makes this Schwarzschild's stellar interior rather than some other static ball.
The parameters are the stellar radius $R$ and $r_s = 2GM/c^2$, the Schwarzschild radius of the whole star, and the areal radius runs over $r \in [0,R]$.

## Step 2. The two candidates

Write

$$f(r) = 1 - \frac{r^2 r_s}{R^3},\qquad f_R = f(R) = 1 - \frac{r_s}{R},$$

so that the two readings of the time term are $g_{tt} = -e^{2\Phi}$ with

$$e^{\Phi}_{\ \text{(line element)}} = \tfrac{1}{2}\left(3\sqrt{f_R} - \sqrt{f}\right),
\qquad
e^{\Phi}_{\ \text{(old inverse)}} = \tfrac{1}{2}\left(3\sqrt{f} - \sqrt{f_R}\right).$$

The second is what the old $g^{tt}$ was the reciprocal of, since $-4/(3\sqrt{f}-\sqrt{f_R})^2$ is exactly $1/g_{tt}$ for that second reading.
The error was therefore a coherent alternative solution and not a slip of a key, so it cannot be dismissed by inspection.

## Step 3. Why the stellar surface does not choose between them

At $r = R$ the two radicals coincide, $\sqrt{f(R)} = \sqrt{f_R}$, and each candidate becomes the same combination of one radical with itself,

$$\tfrac{1}{2}\left(3\sqrt{f_R} - \sqrt{f_R}\right) = \sqrt{f_R},$$

so both give $g_{tt}(R) = -(1 - r_s/R)$ and both join continuously to the exterior Schwarzschild solution.
Moving the factor of three from one radical to the other cannot be seen at a radius where the two radicals are equal.
Both therefore also give the same $g^{tt}(R) = -1/(1-r_s/R)$, which is why the old inverse looked right wherever it was easiest to check.
The surface is silent on the question and the field equations have to be asked instead.

## Step 4. The ansatz and its Einstein tensor

Take the general static spherically symmetric metric in the areal radius,

$$ds^2 = -e^{2\Phi(r)}c^2dt^2 + e^{2\Lambda(r)}dr^2 + r^2\left(d\theta^2 + \sin^2\theta\, d\phi^2\right).$$

Its nonvanishing connection coefficients are

$$\Gamma^t{}_{tr} = \Phi',\qquad
\Gamma^r{}_{tt} = e^{2\Phi-2\Lambda}\Phi',\qquad
\Gamma^r{}_{rr} = \Lambda',$$
$$\Gamma^r{}_{\theta\theta} = -re^{-2\Lambda},\qquad
\Gamma^r{}_{\phi\phi} = -re^{-2\Lambda}\sin^2\theta,$$
$$\Gamma^\theta{}_{r\theta} = \Gamma^\phi{}_{r\phi} = \frac{1}{r},\qquad
\Gamma^\theta{}_{\phi\phi} = -\sin\theta\cos\theta,\qquad
\Gamma^\phi{}_{\theta\phi} = \cot\theta,$$

with the prime a derivative in $r$, and in the $x^0 = ct$ chart $\Gamma^t{}_{tr}$ carries one upper and one lower time index so its factors of $c$ cancel and it is $\Phi'$ there too.
Feeding these through the definitions of Step 1 gives the mixed Einstein tensor

$$G^t{}_t = \frac{e^{-2\Lambda}}{r^2}\left(1 - 2r\Lambda'\right) - \frac{1}{r^2},$$
$$G^r{}_r = \frac{e^{-2\Lambda}}{r^2}\left(1 + 2r\Phi'\right) - \frac{1}{r^2},$$
$$G^\theta{}_\theta = G^\phi{}_\phi = e^{-2\Lambda}\left(\Phi'' + \Phi'^2 - \Lambda'\Phi' + \frac{\Phi' - \Lambda'}{r}\right).$$

## Step 5. The uniform density fixes the radial metric, and then isotropy is linear

The $tt$ equation is the usual mass function, $e^{-2\Lambda} = 1 - 2Gm(r)/(c^2r)$, and a constant $\rho$ makes $m(r) = \tfrac{4}{3}\pi\rho r^3$ grow as $r^3$, so with $2Gm(R)/c^2 = r_s$

$$e^{-2\Lambda} = 1 - \frac{r^2 r_s}{R^3} = f.$$

This is the radial term of the line element and it is not in dispute.
With it, $e^{-2\Lambda}(-2\Lambda'/r) = f'/r = -2r_s/R^3$ and $(f-1)/r^2 = -r_s/R^3$, so

$$G^t{}_t = -\frac{3r_s}{R^3},\qquad
G^r{}_r = \frac{2f\Phi'}{r} - \frac{r_s}{R^3}.$$

The first of these is constant, as a uniform density must make it, and it fixes

$$\rho = \frac{3c^2r_s}{8\pi GR^3} = \frac{3M}{4\pi R^3}.$$

The remaining unknown is $\Phi$, and what determines it is that the fluid is perfect, so the pressure is the same in every direction and $G^\theta{}_\theta = G^r{}_r$.
Using $f\Lambda' = -f'/2 = rr_s/R^3$,

$$G^\theta{}_\theta - G^r{}_r = f\left(\Phi'' + \Phi'^2\right) - \frac{rr_s}{R^3}\Phi' - \frac{f\Phi'}{r},$$

and multiplying by $r$ and using $-r^2r_s/R^3 - f = -1$ this is

$$rf\left(\Phi'' + \Phi'^2\right) - \Phi' = 0.$$

The equation is nonlinear in $\Phi$ but linear in $y = e^{\Phi}$, because $\Phi' = y'/y$ and $\Phi'' + \Phi'^2 = y''/y$:

$$r\left(1 - \frac{r^2r_s}{R^3}\right)y'' = y'.$$

Separating,

$$\frac{y''}{y'} = \frac{R^3}{r\left(R^3 - r^2r_s\right)} = \frac{1}{r} + \frac{rr_s}{R^3 - r^2r_s},$$

so $y' = kr/\sqrt{R^3-r^2r_s}$ and one more integration gives the general solution

$$e^{\Phi} = A + B\sqrt{1 - \frac{r^2 r_s}{R^3}}.$$

Both candidates of Step 2 are of this form, the first with $(A,B) = \left(\tfrac{3}{2}\sqrt{f_R},\,-\tfrac{1}{2}\right)$ and the second with $(A,B) = \left(-\tfrac{1}{2}\sqrt{f_R},\,\tfrac{3}{2}\right)$.
Neither is excluded by the field equations alone, so the question is entirely one of boundary conditions.

## Step 6. The two boundary conditions, and which candidate meets them

Put the general solution into $G^r{}_r$.
With $y = A + B\sqrt{f}$ and $y' = Bf'/(2\sqrt{f}) = -Brr_s/(R^3\sqrt{f})$,

$$\frac{2f}{r}\frac{y'}{y} = -\frac{2Br_s\sqrt{f}}{R^3\left(A + B\sqrt{f}\right)},$$

and therefore

$$\frac{8\pi G}{c^4}p = G^r{}_r = -\frac{r_s}{R^3}\cdot\frac{A + 3B\sqrt{f}}{A + B\sqrt{f}}.$$

There are two conditions to impose, and they are what a ball of fluid of finite radius means.

1. **The star has an edge.** The fluid occupies $r \le R$ and vacuum lies outside it, so the pressure has to reach zero at the surface: $p(R) = 0$, which by that expression for the pressure is $A + 3B\sqrt{f_R} = 0$.
2. **The exterior is Schwarzschild.** The interior joins the vacuum solution of the same mass at $r = R$, so $e^{2\Phi(R)} = 1 - r_s/R$, which is $A + B\sqrt{f_R} = \sqrt{f_R}$.

The first gives $A = -3B\sqrt{f_R}$, and substituting that into the second gives $-2B\sqrt{f_R} = \sqrt{f_R}$, so

$$B = -\tfrac{1}{2},\qquad A = \tfrac{3}{2}\sqrt{f_R},$$

and

$$e^{\Phi} = \frac{1}{2}\left(3\sqrt{1 - \frac{r_s}{R}} - \sqrt{1 - \frac{r^2 r_s}{R^3}}\right).$$

The factor of three belongs to the constant radical.
This is the first candidate of Step 2, the one the line element carries, since the bracket enters $g_{tt}$ squared and $\left(\sqrt{f} - 3\sqrt{f_R}\right)^2 = \left(3\sqrt{f_R} - \sqrt{f}\right)^2$.

The second candidate fails condition 1 and fails it badly.
With $(A,B) = \left(-\tfrac{1}{2}\sqrt{f_R},\,\tfrac{3}{2}\right)$ the numerator at the surface is $A + 3B\sqrt{f_R} = 4\sqrt{f_R}$ while the denominator is $\sqrt{f_R}$, so

$$G^r{}_r\big|_{r=R} = -\frac{4r_s}{R^3},
\qquad
p(R) = -\frac{c^4 r_s}{2\pi GR^3} < 0.$$

It describes a configuration under tension at the surface, held together by something outside itself, which is not a star.
It satisfies condition 2, by Step 3, and that is the whole of its resemblance to the solution.

## Step 7. The pressure and density that come with the answer

With $A = \tfrac{3}{2}\sqrt{f_R}$ and $B = -\tfrac{1}{2}$ the ratio in Step 6 is $3\left(\sqrt{f_R}-\sqrt{f}\right)\big/\left(3\sqrt{f_R}-\sqrt{f}\right)$, and since $8\pi G\rho c^2/c^4 = 3r_s/R^3$,

$$p(r) = \rho c^2\,\frac{\sqrt{1 - \frac{r^2 r_s}{R^3}} - \sqrt{1 - \frac{r_s}{R}}}{3\sqrt{1 - \frac{r_s}{R}} - \sqrt{1 - \frac{r^2 r_s}{R^3}}},
\qquad
\rho = \frac{3M}{4\pi R^3}.$$

This is the textbook interior Schwarzschild pressure.
It is positive throughout $0 \le r < R$, since $f > f_R$ there, it falls monotonically outwards and it vanishes at $r = R$ by construction.

## Step 8. First confirmation: the Newtonian interior potential

Expand in $r_s$, using $\sqrt{1-x} = 1 - x/2 + O(x^2)$, and read off $\Phi = \Phi_{\text{N}}/c^2$ with $r_s = 2GM/c^2$.

For the answer of Step 6,

$$e^{\Phi} \approx 1 - \frac{3r_s}{4R} + \frac{r^2r_s}{4R^3},
\qquad
\Phi_{\text{N}} = -\frac{GM}{2R^3}\left(3R^2 - r^2\right),$$

which is exactly the Newtonian potential inside a uniform sphere, and

$$\nabla^2\Phi_{\text{N}} = \frac{1}{r^2}\frac{d}{dr}\left(r^2\Phi_{\text{N}}'\right) = \frac{3GM}{R^3} = 4\pi G\rho,$$

so it satisfies Poisson's equation with the density of Step 5.

For the rival,

$$e^{\Phi} \approx 1 + \frac{r_s}{4R} - \frac{3r^2r_s}{4R^3},
\qquad
\Phi_{\text{N}} = \frac{GM}{2R^3}\left(R^2 - 3r^2\right),
\qquad
\nabla^2\Phi_{\text{N}} = -\frac{9GM}{R^3}.$$

That is a Newtonian potential sourced by a negative density three times the size of the real one, and it is positive at the centre, so the weak field limit rejects it outright.

## Step 9. Second confirmation: the Buchdahl bound

At the centre the answer of Step 6 gives

$$e^{\Phi(0)} = \frac{1}{2}\left(3\sqrt{1 - \frac{r_s}{R}} - 1\right),$$

which decreases as the star is made more compact and reaches zero at

$$\frac{r_s}{R} = \frac{8}{9}.$$

That is the Buchdahl bound, and the pressure of Step 7 diverges at the centre at exactly the same compactness, since $p(0) = \rho c^2\left(1-\sqrt{f_R}\right)\big/\left(3\sqrt{f_R}-1\right)$ has its pole where $3\sqrt{f_R} = 1$.
Two quantities computed from different components of the same solution failing together at $8/9$ is the known signature of this spacetime.

The rival gives

$$e^{\Phi(0)} = \frac{1}{2}\left(3 - \sqrt{1 - \frac{r_s}{R}}\right),$$

which never vanishes, runs between $1$ and $3/2$, and exceeds one for every star with $r_s > 0$.
A clock at the centre of a star runs slow compared with one at infinity, never fast, so the rival has the central redshift the wrong way round and knows nothing of the Buchdahl bound.

## Step 10. Third confirmation: the stated connection

The stated connection and curvature were already built on the answer of Step 6, which settles the matter internally as well.
From Step 4, $\Gamma^t{}_{tr} = \Phi' = y'/y$, and with $y = \tfrac{1}{2}\left(3\sqrt{f_R} - \sqrt{f}\right)$ and $y' = rr_s/(2R^3\sqrt{f})$,

$$\Gamma^t{}_{tr} = \frac{rr_s}{R^3\sqrt{f}\left(3\sqrt{f_R}-\sqrt{f}\right)} = \frac{rr_s}{R^3\left(3\sqrt{f_Rf} - f\right)}.$$

Now $R^3f = R^3 - r^2r_s$ and $R^3\sqrt{f_Rf} = R\sqrt{(R-r_s)(R^3-r^2r_s)}$, so

$$\Gamma^t{}_{tr} = \frac{rr_s}{-R^3 + r^2r_s + 3R\sqrt{(R-r_s)(R^3-r^2r_s)}},$$

which is the stated value character for character.
The rival does not merely give a different number here, it gives the opposite sign, because $\sqrt{f}$ decreases outwards while $\sqrt{f_R}$ is constant: at $r_s/R = 1/2$ and $r/R = 1/2$ the stated symbol and the answer of Step 6 both give $0.225364609521$ while the rival gives $-0.381958875866$.
Every stated connection coefficient, every curvature component and the Einstein tensor are on the same branch, so the old $g^{tt}$ was the single component that was not.

## Step 11. The corrected component, and $g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu{}_\nu$

The metric is diagonal, so the inverse is componentwise, and with $g_{tt} = -\tfrac{1}{4}\left(3\sqrt{f_R}-\sqrt{f}\right)^2$

$$g^{tt} = -\frac{4}{\left(3\sqrt{1 - \frac{r_s}{R}} - \sqrt{1 - \frac{r^2 r_s}{R^3}}\right)^{\!2}},
\qquad
g^{rr} = 1 - \frac{r^2 r_s}{R^3},
\qquad
g^{\theta\theta} = \frac{1}{r^2},
\qquad
g^{\phi\phi} = \frac{1}{r^2\sin^2\theta}.$$

The last three were already right and are unchanged.
Multiplying the stated inverse metric by the stated metric, both parsed by the checker's own LaTeX reader, gives

$$g^{\mu\alpha}g_{\alpha\nu} = \operatorname{diag}(1,1,1,1) = \delta^\mu{}_\nu$$

identically in $r$, $\theta$, $R$ and $r_s$, with no numerical evaluation and no assumption about the sign of any radical beyond the positivity the parameters already carry.
The $tt$ slot is the only one that moved, and its identity is the algebraic statement

$$-\frac{4}{\left(3\sqrt{f_R}-\sqrt{f}\right)^2}\cdot\left(-\frac{1}{4}\left(3\sqrt{f_R}-\sqrt{f}\right)^2\right) = 1.$$

For contrast, the old component failed this by a factor that depends on position, so it was not a matter of normalisation:

| $r_s/R$ | $r/R$ | old $g^{tt}g_{tt}$ | corrected $g^{tt}g_{tt}$ |
| --- | --- | --- | --- |
| $1/2$ | $1/2$ | $0.319167946686$ | $1$ |
| $4/5$ | $1/4$ | $0.0219507878125$ | $1$ |
| $1/10$ | $9/10$ | $0.959077055123$ | $1$ |
| any | $1$ | $1$ | $1$ |

As in Step 3, the error vanishes at the stellar surface for every mass and every radius, and grows inwards.

## Step 12. What changed

The $tt$ component of the inverse metric changed from

    -\dfrac{4}{\left(3\sqrt{1-\frac{r_s}{R^3}r^2}-\sqrt{1-\frac{r_s}{R}}\right)^2}

to

    -\dfrac{4}{\left(3\sqrt{1 - \frac{r_s}{R}} - \sqrt{1 - \frac{r^2 r_s}{R^3}}\right)^2}

The radicals are also now written in the order and the spacing the line element uses, so the two can be compared by eye.
No other component changed, and in particular the fifty nested radical disagreements the checker reports are left exactly as they were, because they are agreements the simplifier cannot see rather than errors.
Since 22 September 2026 the checker sees them: `norm` splits each radical into the square roots of the irreducible factors of its radicand, which is an identity where every factor is positive, as $R$, $R - r_s$ and $R^3 - r^2 r_s$ are everywhere inside the star.
The fifty are gone, the Einstein tensor and Kretschmann scalar that used to run out of budget now finish, and the whole interior Schwarzschild star checks in about two seconds with no disagreement and nothing unchecked.
