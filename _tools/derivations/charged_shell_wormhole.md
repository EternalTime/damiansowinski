# Schein and Aichelburg's wormhole between charged shells

Sources: F. Schein and P. C. Aichelburg, Phys. Rev. Lett. 77, 4130 (1996), arXiv gr-qc/9606069, read from the arXiv TeX; P. C. Aichelburg and F. Schein, Acta Phys. Pol. B 29, 1025 (1998), read from the journal's PDF.
The equation names below are the labels of the arXiv TeX.

## The charts

All four charts are written by `print_charts.py` through `charged_shell_wormhole.py` beside this file, and each is checked at points before it is written.

- `cartesian`, the founders' (metric), $ds^2 = -c^2dT^2/V^2 + V^2(dx^2 + dy^2 + dz^2)$, with $V$ left free as the live Majumdar-Papapetrou Cartesian chart leaves $U$ free. The check holds it slot by slot to that published chart with $U$ written $V$, and holds its Einstein tensor at two points outside the spheres to the Maxwell stress of $A = c\,dT/V$, in units $G = 4\pi\epsilon_0 = 1$, with $V$ the image series of (6) cut at three images, which is harmonic.
- `bispherical`, the founders' (bisph), with their bispherical coordinates; the founders' focal length $c$ is written $a$, since $c$ is the speed of light. The factor $W = \cosh\mu - \cos\eta$ is held as a name with no rates, so every comparison writes it out, and `bispherical_pretty` writes each value back in $W$, $\sinh\mu$ and the functions of $\eta$, with $\cosh\mu = W + \cos\eta$ and $\sinh^2\mu = (W + \cos\eta)^2 - 1$; the curvature of the coordinates then cancels and the Ricci scalar is one line. The check holds the flat metric pulled back along $x = a\sin\eta\cos\phi/W$, $y = a\sin\eta\sin\phi/W$, $z = a\sinh\mu/W$ to the chart at $V = 1$ at two points, and the bispherical series of $V$ after (bisph) to the Cartesian series (6) at those points, both cut at sixty images, in forty digits.
- `polar`, the founders' (sphcoor), about the centre $z = -d_1$ of one sphere; the check holds the flat metric pulled back along $z = -d_1 + \rho\cos\theta$ to the chart at $V = 1$ at two points.
- `interior`, the founders' (RN) with their charge $e$ written $q$, since $e$ is Euler's number on the page; the check holds it slot by slot to the published spherical chart of `rn_metric` at $r_s = 2m$, $r_q = q$, and its Einstein tensor at three points to a point charge's field, $G^a{}_b = (q^2/r^4)\,\mathrm{diag}(-1, -1, 1, 1)$.

The image charges are the founders' (4) and (5) solved by Lindquist's parameters, $R = a/\sinh\mu_0$, $d_1 = a\coth\mu_0$, $d_n = a\coth n\mu_0$, $m_n = (-1)^{n+1}m_1\sinh\mu_0/\sinh n\mu_0$.
Each shell sits at the areal radius $RV_0 = R + m_1$ inside, the founders' $r_- = RV_0$ with $V_0 = 1 + m_1/R$, and its induced metric is $ds^2 = -c^2d\tau^2 + (R + m_1)^2d\Omega^2$ from both sides.
The shells are in no chart, since the time of one side is not the time of the other.

## What the page does not state

The 1996 letter's numerical critical value, printed as $(R/d_1)_{crit} \approx 9.993858$, cannot be a ratio below one; the 1998 article withdraws it and proves the energy density positive for every $R < d_1$ at small enough $m$ and $|q|$.
The shell pressure is printed with opposite overall signs in the 1996 and 1998 texts, so no pressure formula is given.
