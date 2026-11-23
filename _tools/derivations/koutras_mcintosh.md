# Koutras-McIntosh spacetimes

Conformally flat pure radiation that is not a plane wave: Wils's metric of 1989, Koutras and McIntosh's family of 1996 and Edgar and Ludwig's whole class of 1997.
The charts are written and checked by `koutras_mcintosh.py`, and the drawings by `koutras_mcintosh_drawings.py`, both beside this file.

## The charts and their sources

Every chart is in the signature $(-,+,+,+)$ of Edgar and Ludwig, gr-qc/9612059.

- `koutras_mcintosh`, coordinates $(u, w, x, y)$, is Edgar and Ludwig's eq. (3), Koutras and McIntosh's family with constants $a$ and $b$.
  It is Griffiths and Podolský's line element of their sec. 1 and Kruglikov and Steneker's eq. (1.1), both in the signature $(+,-,-,-)$, with their $v$ for $w$ and, in the second, their $f$ for $2f$.
- `edgar_ludwig`, coordinates $(u, v, x, y)$, is Edgar and Ludwig's eq. (8), with $f$, $g$ and $h$ free functions of $u$; at $g = h = 0$ it is Wils's metric, their eq. (1).
- `kerr_schild`, coordinates $(T, X, Y, Z)$, is Griffiths and Podolský's eq. (12) for Wils's metric with $f = 1/\ell$ constant, turned to $(-,+,+,+)$: the flat metric minus $(x^2 + Y^2)/(2\ell x)$ times the square of $k = (1 + u^2)c\,dT - 2u\,dX - (1 - u^2)dZ$, which is $2x\,du$.
  Its $x$, $u$ and $v$ are those of the Kerr-Schild chart of Kundt's waves, held with the same rates.

Units: $u$ is a pure number and $x$, $y$ lengths; $w$ is a length, $v$ an area, $a$ a pure number, $b$ a length, $f$ an inverse length, $g$ a length and $h$ an area.

## The checks

`check_chart` holds each chart at three points, with test functions of $u$ that are no special case:

- Koutras and McIntosh's chart at $(a, b) = (1, 0)$, $(0, 1)$, $(3/5, 2/7)$ and $(-2/3, 5/4)$ has every component of the Weyl tensor zero and the Ricci tensor $R_{uu} = 4(ax + b)f$ alone, which is Griffiths and Podolský's $\Phi_{22} = 2f/(ax + b)$ in their frame; at $a = 1$, $b = 0$ it is Edgar and Ludwig's eq. (7) with $g = h = 0$.
- Edgar and Ludwig's chart has the Weyl tensor zero and $R_{uu} = 4xf$ alone, and it is their eq. (7) pulled back along $w = v/x$.
- The Kerr-Schild chart is eq. (7) at $f = 1/\ell$, $g = h = 0$ pulled back along Griffiths and Podolský's eq. (7): $w = cT + Z$, $x = \sqrt{X^2 + Z^2 - c^2T^2}$ and $u = (X - x)/(cT + Z)$.

Every curvature scalar of every chart is printed zero, $R = 0$ and $K = 0$, as Koutras and McIntosh state for every scalar built from the Riemann tensor and its derivatives.

## The drawings

- The planes of $u$ and $w$ of Koutras and McIntosh's chart and of $u$ and $v$ of Edgar and Ludwig's draw Wils's metric ($a = 1$, $b = 0$; $g = h = 0$) at $x = \ell$, $y = 0$, under a pulse of radiation switched on and off, $f = e^{-4u^2}/2\ell$: a sandwich wave of the kind Griffiths and Podolský's sec. 5 treats, with $f \ge 0$, so that the energy density is positive.
- The Kerr-Schild chart is drawn as the figure of Kundt's waves, the slice $Y = 0$ of $T$, $X$ and $Z$, since the fronts and their envelope are the same, Griffiths and Podolský's eqs. (8) to (10).
- No conformal diagram is drawn: the amplitude $f(u)$ is free, and the envelope $x = 0$ of the wave fronts is a curvature singularity at which every curvature scalar vanishes (Podolský and Beláň, 2004), so no global completion is established to compactify.
- No embedding diagram is drawn: every wave front is a flat half plane, and no slice is both spacelike and of constant time across the radiation; on the slice $cT = 0$, $Y = 0$ of the Kerr-Schild chart the induced metric is spacelike only where $2r^3 < \ell(r + X)^2$, with $r^2 = X^2 + Z^2$, and the energy density $2f/x$ grows without bound toward the envelope.
- No curvature figure is drawn: both stored scalars, $R$ and $K$, are zero in every chart, so there is nothing for a map to colour.
