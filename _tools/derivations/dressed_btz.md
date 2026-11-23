# The conformally dressed black hole of three dimensions

Cristián Martínez and Jorge Zanelli, "Conformally dressed black hole in 2 + 1 dimensions", Physical Review D 54, 3830 (1996), arXiv gr-qc/9604021, read from its TeX source.

## Step 1. The theory and the solution

The action is their eq. (action), $I = \int d^3x\sqrt{-g}\left[(R + 2l^{-2})/2\kappa - \tfrac12 g^{\mu\nu}\nabla_\mu\Psi\nabla_\nu\Psi - \tfrac{1}{16}R\Psi^2\right]$, with $-l^{-2}$ the cosmological constant.
The field equations are (Eeq), $G_{\mu\nu} - l^{-2}g_{\mu\nu} = \kappa T_{\mu\nu}$, and (Feq), $\Box\Psi - \tfrac18 R\Psi = 0$, with the stress tensor (Tuv).
The stress tensor has no trace, so $R = -6/l^2$, their (R).
Their ansatz (ds) with $N = 1$, $ds^2 = -F\,dt^2 + F^{-1}dr^2 + r^2d\theta^2$, has the solution (Fp), $F = (r + B)^2(r - 2B)/(rl^2) = (r^2 - 3B^2 - 2B^3/r)/l^2$, with the field (psie), $\Psi = \sqrt{8B/(\kappa(r + B))}$, and $B \ge 0$.
The horizon is the null circle $r = 2B$; the Kretschmann scalar is (KR), $12(r^6 + 2B^6)/(l^4r^6)$.
We keep $c$: the time term is $-F\,c^2dt^2$, $t$ a time, and every component is in the chart $x^0 = ct$.

## Step 2. The charts

`static` is (ds) with (Fp), as above.
`eddington_finkelstein_ingoing` is the advanced time chart of their footnote, $ds^2 = -F\,dv^2 + 2\,dv\,dr + r^2d\theta^2$ with $v = t + \int F^{-1}dr$, which with $c$ reads $-F\,c^2dv^2 + 2c\,dv\,dr$, $c\,dv = c\,dt + dr/F$.
Avilés, Maeda and Martínez write the same chart, their eq. (statEF), arXiv 1808.10040.
Martínez's instability paper, gr-qc/9801091 eq. (ds), and Natsuume, Okamura and Sato, hep-th/9910105 eq. (mz), write the static chart with $r_+ = 2B$ in place of $B$; it is the same chart.
The rotating solution of Natsuume and Okamura, hep-th/9911062 (stationaryMZ), is a different spacetime, obtained from this one by a boost that is a coordinate change only locally, and is not published as a chart here.
The minimally coupled frame of Henneaux, Martínez, Troncoso and Zanelli, hep-th/0201170, is a conformally related metric, not a chart of this one, and is not published either.

## Step 3. The checks

`check_chart` in `dressed_btz.py` holds each chart, in forty digits at three points of $(l, B, r)$, outside the horizon, inside it and far from it, to (Eeq) with $\kappa T_{\mu\nu}$ built from $\chi = \sqrt{\kappa}\Psi$, to (Feq), to (R), to (KR), to the two forms of (Fp) and to Martínez's form with $r_+ = 2B$.
At $B = 0$ each chart is the published `btz` chart of the same kind, `stationary` and `eddington_finkelstein_ingoing`, at $M = J = 0$, with $\theta$ for $\phi$ and $l$ for $\ell$; the ingoing chart is also the static one pulled back along $c\,dt = c\,dv - dr/F$.
A planted fault, $r - 3B$ in place of $r - 2B$ in the lapse, fails (Eeq) in both charts, so the check can fail.
