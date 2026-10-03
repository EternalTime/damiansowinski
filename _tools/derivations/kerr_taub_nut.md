# Kerr-Taub-NUT

Kerr's black hole with a NUT parameter, Demiański and Newman's solution of 1966,

$$ds^2 = -\frac{\Delta}{\Sigma}\left(c\,dt - \chi\,d\phi\right)^2 + \frac{\Sigma}{\Delta}dr^2 + \Sigma\,d\theta^2 + \frac{\sin^2\theta}{\Sigma}\left(a\,c\,dt - \left(r^2 + a^2 + l^2\right)d\phi\right)^2,$$

$$\Sigma = r^2 + \left(l + a\cos\theta\right)^2,\qquad \Delta = r^2 - 2mr + a^2 - l^2,\qquad \chi = a\sin^2\theta - 2l\cos\theta.$$

Its five charts are written by `_tools/derivations/print_charts.py --metric kerr_taub_nut`, about six minutes for each chart in $r$ and $\theta$ and ten seconds for Plebański and Demiański's on 2 October 2026, and `verify_metrics.py --system kerr_taub_nut/<chart>` checks each in about three minutes.

## Step 1. Where the metric comes from

The form above is equations (1) and (2) of Bini, Cherubini, Jantzen and Mashhoon (`bini2003`), who take it from Miller (`miller1973`), with their $A$ written $\chi$ and the signature reversed.
It is also equation (1) of Ballon Bordo, Gray, Hennigar and Kubizňák (`bordo2019`) at $s = 0$, and Griffiths and Podolský's Kerr-Newman-NUT metric without charge or acceleration (`griffiths2005`).
The mass enters as the length $m = GM/c^2$, as on the Taub-NUT page, so that $a = 0$ is `taub_nut` term by term, $(c\,dt + 2l\cos\theta\,d\phi)^2$ included, and $l = 0$ is `kerr` with $GM/c^2$ written $m$.
The original paper of Demiański and Newman, Bull. Acad. Polon. Sci. 14, 653 (1966), is indexed by neither Crossref nor INSPIRE, so the page cites what Manko, Martín and Ruiz (`manko2006`), Bini and his collaborators, and Newman himself in his review with Adamo (`adamo2014`) say of it.

## Step 2. The names and how a value is written

Every chart in $r$ and $\theta$ declares $\Sigma$, $\Delta$ and $\chi$ as names for the expressions above.
The checker hands every value back in $r$, $a$, $l$, $m$, $\sin\theta$ and $\cos\theta$, and `KerrTaubNutForms` in `print_charts.py`, a subclass of Kerr-de Sitter's forms, writes it in the names.
It factors the value, writes each factor that is a multiple of $\Sigma$, $\Delta$, $\chi$ or $l + a\cos\theta$ by its name, writes $(r + p)(r - p)$ as $r^2 - p^2$ for $p = l + a\cos\theta$, and writes every sum left over in the shortest of: as it stands in the cosine or the sine, as a polynomial in $r$ and $p$, grouped by the powers of $\chi$, grouped by the powers of $m$, or grouped by the powers of $\Delta$ after $2mr$ is written $r^2 + a^2 - l^2 - \Delta$.
The curvature depends on $r$ and $p$ alone, so every Riemann and Weyl component comes out as a short product, such as $R_{trtr} = -(a^2\sin^2\theta + 2\Delta)(3lr^2p - lp^3 + mr^3 - 3mrp^2)/(\Sigma^3\Delta)$.
Two derivatives of $g_{\phi\phi}$ are written by hand, since no grouping finds them: with $w = r^2 + a^2 + l^2$, $\Sigma + a\chi = w$ and $\partial_\theta\chi = 2p\sin\theta$,

$$\Gamma^r{}_{\phi\phi} = -\frac{\Delta\left(rw\left(\Sigma - a\chi\right)\sin^2\theta - (r - m)\Sigma\chi^2 + r\Delta\chi^2\right)}{\Sigma^3},\qquad \Gamma^\theta{}_{\phi\phi} = -\frac{\left(\Sigma\left(w^2\cos\theta - 2\Delta\chi p\right) + ap\left(w^2\sin^2\theta - \Delta\chi^2\right)\right)\sin\theta}{\Sigma^3}.$$

The factoring names its generators in one order, $r$, $\cos\theta$, $\sin\theta$, $a$, $l$, $m$.
One mixed Riemann component of a Kerr chart once did not factor in forty minutes, which was put down to the order sympy chose; sympy's own order does not depend on the hash seed, and the cause was the random points sympy's factoring lifts at, which `verify_metrics.py` now draws from a fixed seed, `WANG_SEED`.
`by_twist` reduces every power of the cosine above the first by $a\cos^2\theta = a - 2l\cos\theta - \chi$ and never divides polynomials, for the same reason.
`KerrTaubNutForms.by_hand` recognises each of them, and its lowered form, by its value at one rational point, and the chart reads the text back and compares it with the computed value exactly, as it does every printed value.

## Step 3. The curvature

The spacetime is a vacuum, which `kerr_taub_nut_check` holds every chart to before it is written.
It is of Petrov type D with the one Weyl scalar $\Psi_2 \propto (m - il)/(r - i(l + a\cos\theta))^3$ (`bini2003`, equation (4)), so

$$K = \frac{48\left(\left(m^2 - l^2\right)\left(r^6 - 15r^4p^2 + 15r^2p^4 - p^6\right) + 4mlrp\left(3r^4 - 10r^2p^2 + 3p^4\right)\right)}{\Sigma^6},$$

which is $48\,\mathrm{Re}\left((m + il)^2/(r + ip)^6\right)$, Kerr's Kretschmann scalar at $l = 0$ and Taub-NUT's at $a = 0$, and which the script checks against sympy.
It diverges only where $\Sigma = 0$, at $r = 0$ and $\cos\theta = -l/a$, a ring that exists for $|l| \le |a|$ (`griffiths2005`).

## Step 4. The other four charts

Each is checked, slot by slot, to be the first chart pulled back through the Jacobian `kerr_taub_nut_jacobian` writes.

The chart with one Misner string keeps $r$, $\theta$ and $\phi$ and takes $t_N = t + 2l\phi/c$, so that $\chi$ becomes $a\sin^2\theta + 2l(1 - \cos\theta)$, which vanishes on the northern half of the axis, and $r^2 + a^2 + l^2$ becomes $r^2 + (a + l)^2$.
It is Bordo and his collaborators' metric at $s = -n$, Manko and Ruiz's $C = -1$, and the form of Newman, Tamburino and Unti's original metric: the northern half axis is regular and the string lies along $\theta = \pi$.
With $\phi$ periodic in both, the two charts describe spacetimes that agree locally and differ globally, which is Miller's result: for $l \neq 0$ the bifurcation spheres and both halves of the axis cannot be covered in one manifold.

The Kerr charts follow the principal null congruences of that second form, $dv = c\,dt_N + (r^2 + (a + l)^2)\,dr/\Delta$ and $d\tilde\phi = d\phi + a\,dr/\Delta$ for the ingoing one and the opposite signs for the outgoing one.
They are Erbin's (3.12) with his (4.26), $g = (r^2 + a^2 + n^2)/\Delta$ and $h = a/\Delta$ (`erbin2016`), carried to the time $t_N$, which adds $2al$ to $r^2 + a^2 + l^2$.
Then $g_{rr} = 0$ and each chart runs through both roots of $\Delta$ and through $r = 0$, on the half axis as off it.

Plebański and Demiański's chart takes $q = r$, $p = l + a\cos\theta$, $\sigma = -\phi/a$ and $\tau = ct - (a^2 + l^2)\phi/a$, in which

$$ds^2 = -\frac{Q}{p^2 + q^2}\left(d\tau - p^2d\sigma\right)^2 + \frac{p^2 + q^2}{Q}dq^2 + \frac{p^2 + q^2}{P}dp^2 + \frac{P}{p^2 + q^2}\left(d\tau + q^2d\sigma\right)^2,$$

with $Q = q^2 - 2mq + a^2 - l^2$ and $P = a^2 - (p - l)^2$: Griffiths and Podolský's form of Plebański and Demiański's line element without acceleration, charge or cosmological constant, $\epsilon = 1$, $k = a^2 - l^2$ and $n = l$ (`griffiths2005`, `plebanski1976`).
The mass with $q$ and the NUT parameter with $p$ enter alike, $Q = k - 2mq + q^2$ and $P = k + 2lp - p^2$.
$\tau$ is a length and stands as the chart's time, and $\sigma$ is an inverse length.

## Step 5. The values the diagrams are drawn at

$m = 1$, $a = m$ and $l = 5m/4$: Kerr's limiting spin with a twist above it, so that $\Delta = (r - 9/4)(r + 1/4)$, the horizons are $r_+ = 9m/4$ and $r_- = -m/4$, and with $l > a$ no ring singularity is left.
On the equator $g_{tt} = 0$ where $\Delta = a^2$, at $r_E = 1 + \sqrt{41}/4 = 2.601\,m$.
On the regular half axis $f = \Delta/(r^2 + (a + l)^2)$ and $r_* = r + \tfrac{81}{20}\ln|1 - 4r/9| - \tfrac{41}{20}\ln|1 + 4r|$, which `slices.ktn_rstar` sums and `null_rays.py --verify` holds the axis rays to.

The spacetime diagrams draw the principal null rays on the equator in the Boyer-Lindquist chart, in $t$ and $r$ and from above, and in Plebański and Demiański's chart in $\tau$ and $q$, where $d\tau/dq = \pm q^2/Q$; the regular half axis in the chart with one string and in both Kerr charts; and the light cones on the equator in three dimensions.
The Boyer-Lindquist chart draws no plane of $t$ and $r$ on the axis, where it has a Misner string.
The conformal diagram is the tower of the regular half axis, with the chart of one string and both Kerr charts on it.
The equator is not totally geodesic when $l \neq 0$, since $\theta \to \pi - \theta$ is no symmetry, so no conformal diagram of it is drawn; the embedding diagram is its slice of constant Boyer-Lindquist $t$, a surface of the spacetime like any other.
