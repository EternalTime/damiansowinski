# Freund and Rubin's anti-de Sitter space times a 7-sphere

Anti-de Sitter space of four dimensions times a 7-sphere of twice its radius, the solution of supergravity in eleven dimensions with a four-form field on the four large dimensions.
The five charts are written by `_tools/derivations/print_charts.py --metric freund_rubin`, which took 20 seconds for all five on 2 October 2026, and `verify_metrics.py --system freund_rubin/<chart>` checks the five in 16 seconds.

## Step 1. The papers and the source of each chart

Every paper of the candidate entry was verified on Crossref and INSPIRE before anything was built.
Freund and Rubin (1980), Schwarz (1983), Günaydin and Marcus (1985), and Nilsson and Pope (1984) were read in full from the KEK scans of their preprints, EFI 80/35 (`lib-extopc.kek.jp/preprints/PDF/1980/8010/8010222.pdf`), CALT-68-1016 (`.../1983/8308/8308297.pdf`), CALT-68-1185 (`.../1984/8412/8412008.pdf`), and Imperial/TP/83-84/41 (`.../1984/8406/8406336.pdf`); the KEK server answers only a browser's user agent.
Maldacena (1998), Gibbons and Townsend (1993), and Aharony, Gubser, Maldacena, Ooguri, and Oz (2000) were read from arXiv, hep-th/9711200, hep-th/9307049, and hep-th/9905111.
Awada, Duff, and Pope (1983), Kim, Romans, and van Nieuwenhuizen (1985), and Pilch, van Nieuwenhuizen, and Townsend (1984) have no open copy and no KEK scan, and the History says of them only what their abstracts and titles say.
Duff, Nilsson, and Pope's review (Physics Reports 130, 1986) was verified on INSPIRE but could not be read, and is not cited.

Freund and Rubin do not write the metric: their (6) puts the field strength of rank $s$ proportional to the volume form of $s$ dimensions, and their (7a) gives the scalar curvatures of the two factors, $R_{d-s} = (s-1)(d-s)\lambda/(d-2)$ and $R_s = -s(d-s-1)\lambda/(d-2)$.
At $d = 11$ and $s = 4$ each of the seven dimensions carries $\lambda/3$ and each of the four $-2\lambda/3$, the ratio $-2$, so the sphere's radius is twice anti-de Sitter's, which is Maldacena's $R_{sph} = 2R_{AdS}$ of his section 3.2.

- `global`: Aharony et al.'s (2.23), $R^2(-\cosh^2\rho\,d\tau^2 + d\rho^2 + \sinh^2\rho\,d\Omega^2)$, with $\tau = ct/L$, times $4L^2d\Omega_7^2$.
- `conformal`: their (2.24), $\tan\theta = \sinh\rho$, with their $\theta$ written $\chi$.
- `static`: the global chart at $r = L\sinh\rho$.
- `poincare`: their (2.27) with $u = r/L^2$, which is Maldacena's (7.3) with his $U$ written $r$; Maldacena's (3.3) with the 1 of $f$ dropped is the same chart in the radius across the membranes, $r^2 = 4Lu$ of it.
- `proper`: Gibbons and Townsend's (5) at $d = 11$, $p = 2$, $\gamma_x = 2/3$ and $a = 2L$, $e^{2\sigma/L}(-c^2dt^2 + dx^2 + dy^2) + d\sigma^2 + 4L^2d\Omega_7^2$.

The angles of the 7-sphere are $\alpha$, $\beta$, $\gamma$, $\kappa$, $\xi$, $\omega$, and $\psi$; the checker cannot read a Greek letter with a subscript, and the chart printer reads $\zeta$ as sympy's zeta function.

## Step 2. What the script checks

`freund_rubin_check` holds each chart to Freund and Rubin's (1a) with (5a) at $s = 4$, read in the convention $R_{\mu\nu} = R^\alpha{}_{\mu\alpha\nu}$: Weinberg's Einstein equation $R^{\mu\nu} - \frac{1}{2}g^{\mu\nu}R = -8\pi G\theta^{\mu\nu}$ with his opposite sign of the Riemann tensor is $G_{\mu\nu} = 8\pi G\theta_{\mu\nu}$, and in eleven dimensions with $\theta = FF - F^2g/8$ that is $R_{MN} = 8\pi G(F_{MPQR}F_N{}^{PQR} - g_{MN}F^2/12)$.
With $F = f\,\epsilon$ on the anti-de Sitter factor, $F_{\mu PQR}F_\nu{}^{PQR} = -6f^2g_{\mu\nu}$ and $F^2 = -24f^2$, both computed from the components, so the equations hold at $8\pi Gf^2 = 3/4L^2$, with $R^M{}_N = -3/L^2$ on the four large dimensions and $3/2L^2$ on the sphere.
Their (5b) holds since $\sqrt{|g|}F^{0123}$ depends on the angles of the sphere alone, and $F \wedge F$ vanishes as an eight-form on four dimensions.
The Riemann tensor is $-(g_{ac}g_{bd} - g_{ad}g_{bc})/L^2$ on the first block, $(g_{ac}g_{bd} - g_{ad}g_{bc})/4L^2$ on the sphere, and zero across the blocks, component by component.
The static and conformal charts are the global one pulled back along $r = L\sinh\rho$ and $\tan\chi = \sinh\rho$, and the proper distance chart is the Poincaré one along $r = Le^{\sigma/L}$; a difference that `simplify` leaves standing is settled at three points to thirty digits.
The printed Ricci scalar is $-3/2L^2$ and the Kretschmann scalar $117/4L^4$, $24/L^4$ from anti-de Sitter space and $2 \cdot 7 \cdot 6/(2L)^4 = 21/4L^4$ from the sphere.

## Step 3. The drawings

Everything is drawn at $L = 1$.
The spacetime diagrams are the plane of the time and the radial coordinate of each chart, at a point of the 7-sphere, and the plane of $t$ and $\psi$ at $\sigma = 0$, a great circle of the sphere drawn as the arc $2L\psi$; `null_rays.py --verify` holds every ray to $ct \mp r_*$ with $r_* = \arctan\sinh\rho$, $\chi$, $\arctan r$, $-1/r$ and $-e^{-\sigma}$, and to $ct \mp 2\psi$.
A ray goes round the great circle in $4\pi L/c$ and crosses anti-de Sitter space from the boundary to the centre and back out in $\pi L/c$.
The conformal diagram is anti-de Sitter's strip for the global, conformal and static charts, by `ads_global_pq` at $r = \sinh\rho$ and $r = \tan\chi$, and the Poincaré wedge of the strip of AdS2 for the other two, by `poincare_pq(t, 1/r)` with $r = e^\sigma$.
The embedding diagram has three views of $t = 0$: the hyperbolic plane of $\rho$ and $\phi$ in Minkowski space, the cylinder of radius $2L$ of $\sigma$ and $\psi$, and a great 2-sphere of radius $2L$ of $\omega$ and $\psi$.
The spacetime is static, so there is no movie and no stack, and no new kind of diagram.
