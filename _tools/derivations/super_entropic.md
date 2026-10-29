# Super-entropic black holes

Kerr-anti-de Sitter spun up to $a = \ell$, in the chart of Hennigar, Kubizňák and Mann,

$$ds^2 = -\frac{\Delta}{\Sigma}\left(c\,dt - \ell\sin^2\theta\,d\psi\right)^2 + \frac{\Sigma}{\Delta}dr^2 + \frac{\Sigma}{\sin^2\theta}d\theta^2 + \frac{\sin^4\theta}{\Sigma}\left(\ell\,c\,dt - \left(r^2 + \ell^2\right)d\psi\right)^2,$$

$$\Sigma = r^2 + \ell^2\cos^2\theta,\qquad \Delta = \left(\ell + \frac{r^2}{\ell}\right)^2 - 2mr,\qquad \psi \sim \psi + \mu.$$

Its two charts are written by `_tools/derivations/print_charts.py --metric super_entropic`, which took 24 seconds on 29 October 2026, and `verify_metrics.py --system super_entropic/<chart>` checks each in seconds.

## Step 1. Where the metric comes from

The line element is equations (8) to (10) of Hennigar, Kubizňák and Mann (`hennigar2015prl`) and (2.3) to (2.4) of Hennigar, Kubizňák, Mann and Musoke (`hennigar2015jhep`), with the charge $q = 0$ and their $l$ written $\ell$.
They reach it from Carter's chart of Kerr-Newman-anti-de Sitter in the form of Hawking, Hunter and Taylor-Robinson, the chart `boyer_lindquist` of `kerr_de_sitter` with $\Lambda = -3/\ell^2$ and $r_s = 2m$, by writing $\psi = \phi/\Xi$ with $\Xi = 1 - a^2/\ell^2$ and letting $a \to \ell$.
Every term then has a finite limit: $a\sin^2\theta\,d\phi/\Xi = a\sin^2\theta\,d\psi \to \ell\sin^2\theta\,d\psi$, $(r^2 + a^2)\,d\phi/\Xi \to (r^2 + \ell^2)\,d\psi$, $\Delta_\theta = 1 - (a^2/\ell^2)\cos^2\theta \to \sin^2\theta$, and $\Delta_r = (r^2 + a^2)(1 + r^2/\ell^2) - 2mr \to (r^2 + \ell^2)^2/\ell^2 - 2mr$, which is $\Delta$.
The period of $\phi$, $2\pi$, would make the period of $\psi$ equal to $2\pi/\Xi$, which diverges, so the limit leaves $\psi$ on the whole real line and the authors close it again with a period $\mu$ of their choosing.

## Step 2. The Carter-Plebański chart

Klemm (`klemm2014`, his (1) and (6)) writes the same geometry as Carter and Plebański's metric with the quartic $P(p)$ given two double roots, $P = (p^2 - \ell^2)^2/\ell^2$, and zero charges, and shifts $\tau \to \tau + \ell^2\sigma$ so that the boundary has no closed timelike curves with $\sigma$ periodic.
With the shift made his line element is

$$ds^2 = -\frac{Q}{q^2 + p^2}\left(c\,d\tau + \left(\ell^2 - p^2\right)d\sigma\right)^2 + \frac{q^2 + p^2}{Q}dq^2 + \frac{q^2 + p^2}{P}dp^2 + \frac{P}{q^2 + p^2}\left(c\,d\tau + \left(\ell^2 + q^2\right)d\sigma\right)^2,$$

with $Q = (\ell + q^2/\ell)^2 - 2mq$.
Hennigar, Kubizňák and Mann's (11) relates the two: $\tau = t$, $q = r$, $p = \ell\cos\theta$, $\sigma = -\psi/\ell$, and $\sigma$'s period is $\mu/\ell$, Klemm's $L$.
Then $c\,dt - \ell\sin^2\theta\,d\psi = c\,d\tau + (\ell^2 - p^2)\,d\sigma$, $\ell\,c\,dt - (r^2 + \ell^2)\,d\psi = \ell(c\,d\tau + (\ell^2 + q^2)\,d\sigma)$, $P = \ell^2\sin^4\theta$, and $dp^2/P = d\theta^2/\sin^2\theta$, so the two line elements agree term by term.
`super_entropic_check` holds the first chart, slot by slot, to the second pulled back through that change, whose Jacobian is diagonal, $\mathrm{diag}(1, 1, -\ell\sin\theta, -1/\ell)$.

## Step 3. The names and the curvature

The first chart declares $\Sigma$ and $\Delta$ and the second $Q$ and $P$ as names for the expressions above.
`SuperEntropicForms` in `print_charts.py` is Kerr-de Sitter's factoring with these names: a factor that is a multiple of $\Sigma$, of $\ell^2\Delta$ or of $\ell^2Q$ is written by its name, and a sum left over in the shortest of its forms as it stands, in the sine, or grouped by the powers of $m$.
The spacetime is an Einstein space, $R_{\mu\nu} = -(3/\ell^2)g_{\mu\nu}$ and $R = -12/\ell^2$, which `super_entropic_check` holds both charts to before either is written.
The curvature scalars do not see the change $\psi = \phi/\Xi$, and Kerr-anti-de Sitter's are continuous in $a$, so the Kretschmann scalar is Kerr-de Sitter's at $a = \ell$, $r_s = 2m$ and $\Lambda = -3/\ell^2$:

$$K = \frac{48m^2\left(r^2 - \ell^2\cos^2\theta\right)\left(\Sigma^2 - 16\ell^2r^2\cos^2\theta\right)}{\Sigma^6} + \frac{24}{\ell^4},$$

which the script checks against sympy in both charts.
It diverges only where $\Sigma = 0$, at $r = 0$ on the equator, which lies inside the inner horizon and outside both charts.

## Step 4. Domains

Both charts cover the region outside the event horizon $r_+$, the largest root of $\Delta$.
The poles $\theta = 0, \pi$, where $p = \pm\ell$, are left out: Hennigar, Kubizňák, Mann and Musoke show that no null geodesic from the bulk reaches them in a finite affine parameter and that they lie at infinite proper distance along every slice of constant $t$ and $r$ (their section II.B, footnote 1, and II.C), so $\theta$ runs over the open interval $(0, \pi)$.
Horizons exist for $m \ge m_0 = 2r_0(r_0^2/\ell^2 + 1)$ with $r_0 = \ell/\sqrt{3}$, their (2.6) and (2.7) at $q = 0$, which is $m_0 = 8\ell/(3\sqrt{3})$; at $m = m_0$ the two roots meet and the hole is extremal.
