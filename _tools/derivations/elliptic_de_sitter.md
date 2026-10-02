# Elliptic de Sitter space

The five charts of `elliptic_de_sitter.json` are written by `print_charts.py --metric elliptic_de_sitter`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, the identification in each, and what the diagrams draw.

## Step 1. The spacetime

De Sitter space is the hyperboloid $-X_0^2 + X_1^2 + X_2^2 + X_3^2 + X_4^2 = \ell^2$ in flat space of five dimensions, with $\Lambda = 3/\ell^2$.
Elliptic de Sitter space is that hyperboloid with $X$ and $-X$ one event (`schrodinger1956`, chapter I section 3; `parikh2003`, section 3).
The map $X \to -X$ has no fixed point on the hyperboloid and commutes with every Lorentz transformation of the five dimensions, so the quotient is smooth and keeps all ten Killing vectors.
Every local quantity is de Sitter's: $R_{\mu\nu} = (3/\ell^2)g_{\mu\nu}$, no Weyl tensor, $R = 12/\ell^2$ and $K = 24/\ell^4$.

The charts are the five of `spradlin2003`, section 2, about the observer at the pole $X_1 = X_2 = X_3 = 0$, $X_4 > 0$, with $n$ the unit vector of $\theta$ and $\phi$:

- `global`: $X_0 = \ell\sinh(ct/\ell)$, $X_4 = \ell\cosh(ct/\ell)\cos\chi$, $X_i = \ell\cosh(ct/\ell)\sin\chi\,n_i$.
- `conformal`: $\tan\eta = \sinh(ct/\ell)$, so $X_0 = \ell\tan\eta$, $X_4 = \ell\cos\chi/\cos\eta$, $X_i = \ell\sin\chi\,n_i/\cos\eta$.
- `kruskal`: $X_0 = \ell(U + V)/(1 - UV)$, $X_4 = \ell(V - U)/(1 - UV)$, $X_i = \ell(1 + UV)\,n_i/(1 - UV)$, the answer to the lectures' Exercise 1 with the observer's static patch at $U < 0 < V$.
- `static`: $X_0 = \sqrt{\ell^2 - r^2}\sinh(ct/\ell)$, $X_4 = \sqrt{\ell^2 - r^2}\cosh(ct/\ell)$, $X_i = r\,n_i$.
- `planar`: $X_0 = \ell\sinh(ct/\ell) + e^{ct/\ell}\rho^2/2\ell$, $X_4 = \ell\cosh(ct/\ell) - e^{ct/\ell}\rho^2/2\ell$, $X_i = e^{ct/\ell}x_i$, with $\rho^2 = x^2 + y^2 + z^2$, the chart of the expanding half $X_0 + X_4 > 0$.

`elliptic_de_sitter_check` holds each chart to the hyperboloid's metric pulled back along these, to $R_{\mu\nu} = (3/\ell^2)g_{\mu\nu}$ and to a vanishing Weyl tensor.

## Step 2. The identification in each chart

In the global chart $X \to -X$ is $t \to -t$, $\chi \to \pi - \chi$, $\theta \to \pi - \theta$, $\phi \to \phi + \pi$ (`parikh2003`, section 3), in the conformal chart the same with $\eta \to -\eta$, and in the Kruskal chart $U \to -U$, $V \to -V$ with the same turn of the sphere.
`EDS_ANTIPODE` in `print_charts.py` holds the three maps, and the check holds each to changing the sign of every $X_I$.
A fundamental domain is half of space for all time, $\chi \le \pi/2$, or $V \ge U$, which is the same half since $\tan\chi = (1 + UV)/(V - U)$.
Its edge $\chi = \pi/2$ is carried onto itself with the time reversed, which each chart states as its last domain.

The static patch is $X_4 > |X_0|$, since $X_4^2 - X_0^2 = \ell^2 - r^2$, and the planar chart is $X_0 + X_4 = \ell e^{ct/\ell} > 0$.
The antipodal map carries each off itself, so neither holds an antipodal pair and neither states an identification.
The antipode of the static patch is the static patch of the opposite observer, so in the quotient the two are one patch, with $t$ reversed if both times are taken to increase with $X_0$.
The planar chart and its antipode are the two halves of the hyperboloid, so in the quotient the planar chart covers everything but the null surface $X_0 + X_4 = 0$.

De Sitter's own proposal of 1917 is this identification (`desitter1917`, article 2).
His static chart with $r = \ell\sin\chi$ has $X_0 = \ell\cos\chi\sinh(ct/\ell)$ and $X_4 = \ell\cos\chi\cosh(ct/\ell)$, and opposite points of the sphere of space at one $t$, $\chi \to \pi - \chi$ with $n \to -n$, are $X$ and $-X$.
`_tools/test_elliptic_de_sitter.py` holds that.

## Step 3. What the identification does

With $Z = X\cdot Y/\ell^2$, two events are joined by a geodesic of the hyperboloid when $Z > -1$, and by none when $Z < -1$: the geodesics through $X$ are $X\cos s + V\sin s$, $X\cosh s + V\sinh s$ and $X + sV$, on which $Z$ is $\cos s$, $\cosh s$ and $1$.

- $Z(X, -X) = -1$ and $(X - (-X))^2 = 4\ell^2$: antipodes are spacelike separated.
- An event on the light cone of $X$ has $Z = 1$ and one on the light cone of $-X$ has $Z = -1$, so the two cones share no event and the quotient has no closed timelike curve (`parikh2003`, section 3.1).
- $Z(X, -Y) = -Z(X, Y)$, so where no geodesic runs from $X$ to $Y$ a timelike one runs to $-Y$: any two events of the quotient are joined, as Calabi and Markus proved of every quotient that is not time orientable (`calabi1962`, as `hawking1973` section 5.2 states it).
- $X \to -X$ carries a tangent vector $V$ to $-V$, and a timelike vector is future directed when $V_0 > 0$, so the map reverses the direction of time and the quotient has none.
- The past of the whole world line of the observer at the pole is $X_4 > X_0$, which holds for exactly one of $X$ and $-X$ off the null surface $X_4 = X_0$: the observer sees one event of every antipodal pair (`schrodinger1956`, the theorem of the footnote to section 3).

The test file holds each of these at random events.

## Step 4. The diagrams

- Spacetime diagrams: the global and conformal charts along the line through the observer, $\chi$ from $-\pi/2$ to $\pi/2$, whose two edges are one line of events, the right edge at $t$ being the left edge at $-t$, a Möbius band; the Kruskal chart on $V \ge U$, with both horizons marked; the static chart on its plane of $t$ and $r$; and the planar chart on its plane of $t$ and $x$. `--verify` checks the rays against $\mathrm{gd}(ct/\ell) \pm \chi$ with $\mathrm{gd}(s) = 2\arctan\tanh(s/2)$, $\eta \pm \chi$, $U$ and $V$, $ct \pm \ell\,\mathrm{artanh}(r/\ell)$, and $x \mp \ell e^{-ct/\ell}$.
- Conformal diagram: half of de Sitter's square, $0 \le \chi \le \pi/2$, in five views. The Kruskal chart enters as $p = \arctan V - \pi/4$, $q = \arctan U + \pi/4$. Each view carries the glued edge, one event $P$ of it drawn at both of its places, and a light ray from an event $E$ beyond the observer's future horizon out through $P$ and in to the observer. The planar view tints the part carried in by the identification apart from the part that lies in the half square directly.
- Embedding diagram: the moment $t = 0$ of the global chart on the equatorial plane, the hemisphere $\rho = \ell\sin\chi$, $z = -\ell\cos\chi$, whose rim is glued to itself point to opposite point, a real projective plane. Two opposite meridians are marked as one geodesic through the observer, closed after the length $\pi\ell$, de Sitter's "total length" of a straight line in elliptical space.

## Sources as read

Schrödinger's book is lending only at archive.org; its preface and pages 7 to 14 were read through the full text search of the scan, phrase by phrase, and every quotation in the history was matched word for word.
Hawking and Ellis, pages 130 and 131, were read the same way.
Crossref gives N. Sánchez as the one author of `sanchez1987`, as her own later papers cite it; INSPIRE's record adds A. Folacci.
