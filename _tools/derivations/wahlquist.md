# Wahlquist's rotating fluid

Why each chart of `wahlquist` is the one published, where it comes from, and what `print_charts.py` holds it to.
The builder is `wahlquist` in `print_charts.py`, the forms are `WahlquistForms`, and the checks are `wahlquist_check`.

## Step 1: the sources

Wahlquist's paper is H. D. Wahlquist, *Physical Review* **172**, 1291 (1968), doi:10.1103/PhysRev.172.1291, which is not open.
Its line element is taken as two open papers reproduce it from that source: Cuchí, Martín, Molina and Ruiz, arXiv:1301.4962, section 3, in Wahlquist's own letters $\xi$, $\eta$, $h_1$, $h_2$, $b$, $k$, $r_0$, $c$ and $\eta_0$, and Bradley, Fodor, Marklund and Perjés, arXiv:gr-qc/9910001, section 3, who write $\zeta$, $\xi$ and $\kappa$ for his $\xi$, $\eta$ and $b$.
The two agree, and Stephani, Kramer, MacCallum, Hoenselaers and Herlt give the same form as their (21.57).
The check below then holds that line element to the perfect fluid it is said to be, so the reproduction is tested and not only copied.

Mars's two charts are M. Mars, *Physical Review D* **63**, 064022 (2001), arXiv:gr-qc/0101021: his (9), the form in which every function is regular at $\beta = 0$, and his (22), the extension across $V = 0$.
Hinoue, Houri, Rugina and Yasui, arXiv:1402.6904, write the first as their (2.1).
Whittaker's sphere is J. M. Whittaker, *Proceedings of the Royal Society A* **306**, 1 (1968), in the form Bradley and his coauthors give it in their section 3, with $\sin X = kr/r_0$.

## Step 2: Wahlquist's chart

$$ds^2 = -\frac{h_1 - h_2}{\xi^2 + \eta^2}\left(c\,dt - \gamma r_0\left(\frac{\xi^2 h_2 + \eta^2 h_1}{h_1 - h_2} - \eta_0^2\right)d\phi\right)^2 + r_0^2\left(\xi^2 + \eta^2\right)\left(\frac{d\xi^2}{(1 - k^2\xi^2)h_1} + \frac{d\eta^2}{(1 + k^2\eta^2)h_2} + \frac{\gamma^2 h_1 h_2}{h_1 - h_2}d\phi^2\right)$$

The signature is the collection's, $(-,+,+,+)$, where Wahlquist's is the opposite.
His constant $c$ is written $\gamma$, since $c$ is the speed of light on every page.
It multiplies $d\phi$, so that $\phi$ is periodic in $2\pi$ with no cone on the axis: near $\eta = \eta_0$ the proper distance from the axis is $2r_0\sqrt{(\xi^2 + \eta_0^2)(\eta_0 - \eta)/((1 + k^2\eta_0^2)|h_2'|)}$ and the circle of $\phi$ has the radius $\gamma r_0\sqrt{(\xi^2 + \eta_0^2)|h_2'|(\eta_0 - \eta)}$, which agree when $2/\gamma = \sqrt{1 + k^2\eta_0^2}\,|h_2'(\eta_0)|$.
The constant $\eta_0^2$ in the dragging term makes $g_{t\phi}$ vanish on the axis.
Both $\eta_0$ and $\gamma$ are parameters of the chart, since neither has a closed form, and the local geometry is the fluid's for any value of either.

The chart is the oblate spheroidal chart of flat space when the fluid is taken away: $\xi = 0$ is a disc, $\eta = 0$ the equatorial plane outside it, and $\xi = \eta = 0$ the ring where they meet, on which $(h_1 - h_2)/(\xi^2 + \eta^2)$ is $0/0$ as written and $1$ in the limit.
The centre of the body is $\xi = 0$, $\eta = \eta_0$.

## Step 3: Mars's chart, and why two functions are held

$$ds^2 = -\frac{V}{v_1 + v_2}\left(c\,d\tau - v_1 d\sigma\right)^2 + \frac{U}{v_1 + v_2}\left(c\,d\tau + v_2 d\sigma\right)^2 + \left(v_1 + v_2\right)\left(\frac{dy^2}{V} + \frac{dz^2}{U}\right)$$

with $v_1 = \sinh^2(\beta z)/\beta^2$ and $v_2 = \sin^2(\beta y)/\beta^2$.
Mars's $U$ and $V$ are written here through a NUT function $n = a_1 + \mu_0 z/\beta^2$ and a mass function $m = a_2 - \mu_0 y/\beta^2$, as $U = Q_0 - (\nu_0 + \mu_0/\beta^2)v_1 + n\sinh(2\beta z)/2\beta$ and $V = Q_0 + (\nu_0 + \mu_0/\beta^2)v_2 + m\sin(2\beta y)/2\beta$, which is his (8) regrouped.
Each holds its coordinate bare beside a sine of it, so its derivatives are algebraic in the function itself: $\partial_z^2 U = 4\beta^2(U - Q_0) - 2\nu_0 + 4\mu_0 v_1$, and $\sinh(2\beta z)\,\partial_z U$ is a polynomial in $U$ and $v_1$.
`HELD` in `verify_metrics.py` therefore holds $U$ and $V$ as functions while the tensors are built, as it holds Wahlquist's $h_1$ and $h_2$ and Whittaker's $F$, and `WahlquistForms.reduce` writes each derivative and each bare coordinate in the function.
Built with $U$ and $V$ written out, the Einstein tensor alone took 112 seconds on 2 October 2026; held, the whole chart prints in 108.

## Step 4: the Weyl tensor is two scalars

The solution is of Petrov type D, with repeated principal null directions $e_0 \pm e_1$ in the frame of the line element, $e^0 \propto c\,d\tau - v_1 d\sigma$, $e^1 \propto dy$, $e^2 \propto dz$ and $e^3 \propto c\,d\tau + v_2 d\sigma$.
Its Weyl tensor then has two frame components, $W_1 = C_{\hat 0\hat 1\hat 0\hat 1}$ and $W_2 = C_{\hat 0\hat 1\hat 2\hat 3}$, and every other is one of them times $\pm 1$ or $\pm\tfrac{1}{2}$.
In $n$ and $m$ the constants $Q_0$ and $\nu_0$ drop out of both:

$$W_1 = -\frac{n\sinh(2\beta z)A_2 + m\sin(2\beta y)A_1}{2\beta(v_1 + v_2)^3} - \frac{\mu_0\left(\beta^2(v_1^2 - 4v_1v_2 + v_2^2) + 3v_1 - 3v_2\right)}{3\beta^2(v_1 + v_2)^2}$$

$$W_2 = \frac{n\sin(2\beta y)A_1 - m\sinh(2\beta z)A_2}{2\beta(v_1 + v_2)^3} - \frac{\mu_0\sinh(2\beta z)\sin(2\beta y)}{2\beta^4(v_1 + v_2)^2}$$

with $A_1 = 2\beta^2v_1(v_1 - v_2) + 3v_1 - v_2$ and $A_2 = 2\beta^2v_2(v_1 - v_2) - v_1 + 3v_2$.
At $\beta = 0$ and $\mu_0 = 0$ they are the real and imaginary parts of the Kerr-NUT $\Psi_2$.
Each chart names $W_1$ and $W_2$ among its parameters, and `WahlquistForms.prepare` writes every Weyl component as $W_1$ and $W_2$ times the products of the frame, and every Riemann component as that plus the part the Ricci tensor fixes.
No step trusts these formulas: each printed component is read back by the checker's reader, with the names written out, and compared with the tensor built from the line element.
The Kretschmann scalar is $12(W_1^2 - W_2^2)$ plus the fluid's part, checked the same way.
The connection is printed with $\partial_z U$ and $\partial_y V$ standing, which is how Carter's separable metrics are short.

## Step 5: the root of $1 - k^2\xi^2$

The checker's canonical form splits a radical into the roots of the irreducible factors of its radicand, and writes $\sqrt{1 - k^2\xi^2}$ as $i\sqrt{k\xi - 1}\sqrt{k\xi + 1}$.
That is the same algebra on the other branch of the root: every identity among the tensors holds on both, so every comparison the checker makes is sound, and a number taken from a canonical form is the wrong one.
So `WahlquistForms.polynomial` reads the canonical roots back as the one radical the line element writes, the printer writes that radical as one symbol, and every numerical check, the pullback from Mars's chart and the limit to Whittaker's, takes $h_1$ and $h_2$ as the reader reads them and never from a canonical form.
The diagrams read the file through the reader and are not affected.

## Step 6: a derivative of a function with a subscript

The reader took a `\partial` only of a function whose name carries no subscript.
`expand_partials` now takes the numeral subscript with the name where the two together are a declared function, so that `\partial_\xi h_1` is the derivative of $h_1$; `_tools/test_reader_names.py` holds both readings.
`null_rays.load` holds every function of `HELD`, so the drawings read such a derivative too.

## Step 7: what each chart is held to

- Every chart: $G^\mu{}_\nu$ is a perfect fluid's moving along the Killing vector of the time, $\rho + 3p$ is the stated constant, and the pressure is Wahlquist's $p = \tfrac{1}{2}\mu_0(1 - b^2 f)$, in Mars's constants $p = \mu_0 + \beta^2 g_{\tau\tau}$ with $8\pi G/c^4 = 1$.
- Wahlquist's chart is Mars's pulled back along $y = r_0\arcsin(k\xi)/k$, $z = r_0\,\mathrm{arsinh}(k\eta)/k$, $c\tau = ct + \gamma r_0\eta_0^2\phi$ and $\sigma = \gamma\phi/r_0$, at $Q_0 = r_0^2$, $\nu_0 = 1$, $a_1 = a_2 = 0$, $\beta = k/r_0$ and $\mu_0 = k^2/(br_0)^2$, at three points in forty digits. Mars's $\mu_0$ is half of Wahlquist's, as a curvature.
- Mars's ingoing chart is his canonical one pulled back along $d\tau = dv - v_2\,dy/V$ and $d\sigma = -d\phi + dy/V$, exactly.
- Whittaker's chart is the limit $k \to 0$ of Wahlquist's at $r_0 = kR_0$, $\xi = \sin X/k$ and $\eta = \cos\theta$, taken at $k = 10^{-20}$ in sixty digits and held to $10^{-15}$, since $g_{t\phi}$ vanishes with the first power of $k$.

## Step 8: the body the diagrams draw

$k = 3/10$ and $b = 4/5$ in units of $r_0$: $\eta_0 = 1.02460115214$, $\gamma = 1.02881004732$, and the surface of zero pressure crosses the equatorial plane at $\xi = 2.80977761761$, where $k\xi = \sin X_s$ with $X_s\cot X_s = b^2$, the surface of Whittaker's sphere, and the axis at $\xi = 2.9124$.
At $k = 1/2$ and above with this $b$, $h_2$ has no root and there is no axis.
From the centre the pole is $3.332\,r_0$ along the axis and the equator $3.310\,r_0$ across the disc and the plane, so this body is prolate; `_tools/test_wahlquist.py` holds those numbers.
No conformal diagram is drawn: the fluid ends on its surface and no exterior is known that carries it to infinity.
