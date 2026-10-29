# Buchdahl's fluid sphere: the chart, the source and what the diagrams draw

Buchdahl's sphere is a static star of perfect fluid whose density falls from the centre and reaches zero only at infinity, the relativistic counterpart of the Newtonian polytrope of index 5.
This note records where the chart of `buchdahl_sphere.json` comes from, how the checker reads its held name, and what the diagrams take from it.

## Step 1. The equations of the paper

Buchdahl (1964) writes isotropic coordinates with signature $(+,-,-,-)$, his (2.1),

$$ds^2 = -e^{\mu(r)}\left(dr^2 + r^2d\Omega^2\right) + c^2e^{\nu(r)}dt^2,$$

and sets, his (2.2),

$$e^\mu = (1 + f)^4, \qquad e^\nu = (1 - f)^2(1 + f)^{-2},$$

the form Schwarzschild's exterior takes with $f = GM/(2c^2r)$, his (2.3).
Equating the two expressions for the pressure, his (2.4) and (2.5), leaves his (2.7), $ff'' - 3f'^2 - r^{-1}ff' = 0$, solved by his (2.8),

$$f = \frac{a}{2c^2}\left(1 + kr^2\right)^{-1/2},$$

with the pressure and density of his (2.9) and (2.10), $p = bc^2f^6(1 + f)^{-5}(1 - f)^{-1}$, $\rho = 3bf^5(1 + f)^{-5}$, $b = (8/\pi)G^{-1}a^{-4}kc^{10}$.
Far out $f \to a/(2c^2\sqrt{k}\,r)$, so the mass is his (2.14), $M = a/(G\sqrt{k})$.
The paper was read in full from the scan the Astrophysics Data System holds.

## Step 2. The chart

With $r_s = 2GM/c^2$ and $R = k^{-1/2}$, Igata, Omamiuda and Takamori's core radius (their (14)), $a = c^2r_s/(2R)$ and

$$f = \frac{r_s}{4\sqrt{r^2 + R^2}},$$

which is Schwarzschild's $r_s/(4r)$ at $R = 0$.
The file turns the signature around and writes

$$ds^2 = -\frac{(1 - f)^2}{(1 + f)^2}c^2dt^2 + (1 + f)^4\left(dr^2 + r^2d\Omega^2\right).$$

The centre has $f_c = r_s/(4R)$, and $g_{tt}$ vanishes nowhere while $f_c < 1$, that is $r_s < 4R$.
Buchdahl's bound $q_c = 3p_c/(\rho_cc^2) \le 1$, his (2.15), is $f_c \le 1/2$, $r_s \le 2R$, by his (2.12) $q = f/(1 - f)$.
The dominant energy condition holds everywhere for $f_c \le 3/4$, $r_s \le 3R$ (Igata, Omamiuda and Takamori, their (8) and (9)).
The areal radius is $r(1 + f)^2$.

## Step 3. How the checker reads $f$

$f$ is held as a function of $r$ (`HELD` in `verify_metrics.py`) with its slope declared (`RATES`), $f' = -rf/(r^2 + R^2)$.
The slope fixes $f$ up to a constant factor, which is $r_s$, and the field equations hold for every value of it, so a comparison made between rational functions of $f$, $r$ and $R$ can disagree falsely and never agree falsely.
Every printed value is such a rational function, and $r_s$ stands only in the definition of $f$.

## Step 4. The check

`buchdahl_sphere_check` in `print_charts.py` compares, at random rational points and to forty digits:

- the declared slope with the derivative of the definition of $f$;
- $G^t{}_t$ with $-8\pi G\rho/c^2 = -3072R^2f^5/(r_s^4(1 + f)^5)$, Buchdahl's (2.9) and (2.10) in these names;
- $G^r{}_r$, $G^\theta{}_\theta$ and $G^\phi{}_\phi$ with $8\pi Gp/c^4 = 1024R^2f^6/(r_s^4(1 - f)(1 + f)^5)$;
- $3p/(\rho c^2)$ with $f/(1 - f)$, his (2.12);
- the Ricci tensor at $R = 0$ with zero, Schwarzschild's vacuum.

Printed, $8\pi G\rho/c^2 = 12fR^2/((r^2 + R^2)^2(1 + f)^5)$ and $8\pi Gp/c^4 = 4f^2R^2/((r^2 + R^2)^2(1 - f)(1 + f)^5)$, the same values with $(r^2 + R^2)^2 = r_s^4/(256f^4)$.
