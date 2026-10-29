# The eye of the storm

Simpson and Visser's eye of the storm, which is Ghosh's nonsingular rotating black hole, is Kerr's metric with the mass $m$ replaced by $me^{-\ell/r}$ wherever it appears.
In the coordinates of Boyer and Lindquist and in the ingoing Kerr chart it is

$$ds^2 = -\frac{\Delta}{\Sigma}\left(c\,dt - a\sin^2\theta\,d\phi\right)^2 + \frac{\Sigma}{\Delta}dr^2 + \Sigma\,d\theta^2 + \frac{\sin^2\theta}{\Sigma}\left(\left(r^2 + a^2\right)d\phi - a\,c\,dt\right)^2,$$

$$ds^2 = -\frac{\Delta}{\Sigma}\left(dv - a\sin^2\theta\,d\tilde\phi\right)^2 + 2\,dr\left(dv - a\sin^2\theta\,d\tilde\phi\right) + \Sigma\,d\theta^2 + \frac{\sin^2\theta}{\Sigma}\left(\left(r^2 + a^2\right)d\tilde\phi - a\,dv\right)^2,$$

with $\Sigma = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 + a^2 - 2mr\,e^{-\ell/r}$.

Both charts are written by `_tools/derivations/print_charts.py --metric eye_of_the_storm`, and `verify_metrics.py --system eye_of_the_storm/boyer_lindquist` and `--system eye_of_the_storm/kerr_ingoing` check each.
The sources are A. Simpson and M. Visser, JCAP 2022(03) 011, arXiv:2111.12329 (v2 read in full), and S. G. Ghosh, Eur. Phys. J. C 75, 532 (2015), arXiv:1408.5668 (v4 read in full).
The literature review behind the entry, with every claim, query and quotation, is kept with the task's records, not in the repository.

## Step 1. The metric and its two charts

Simpson and Visser's (2.26)-(2.27) is the Boyer-Lindquist line element above with $c = 1$, and Ghosh's (1) is the same metric multiplied out, with his $M$ and $k$ our $m$ and $\ell$.
Ghosh's (13)-(14) carry it to the ingoing Kerr chart by $dv = c\,dt + (r^2 + a^2)\,dr/\Delta$ and $d\tilde\phi = d\phi + a\,dr/\Delta$; with $M$ and $k$ constant his (14) is the second line element above.
`eye_of_the_storm_check` in `print_charts.py` holds the ingoing chart to Ghosh's (14) component by component and to the Boyer-Lindquist chart pulled back along that Jacobian, before anything is written.
Every principal null ray of the ingoing family is a line of constant $v$, $\theta$ and $\tilde\phi$ in it, as in Kerr's, since the substitution leaves the principal one form $c\,dt - a\sin^2\theta\,d\phi$ and its pairing with $dr$ unchanged.

## Step 2. Units

$\ell/r$ is a number, so $\ell$ is a length, and $m = GM/c^2$ and $a = J/(Mc)$ are lengths as in Kerr's entry.
`DIMENSIONS` declares $[m] = [a] = [\ell] = L$ and $[\Sigma] = [\Delta] = L^2$, and the dimensional pass balances every published term.
The reader reads $e^{-\ell/r}$ as the exponential, and it is written $2mr\,e^{-\ell/r}$ with a thin space, since $mre$ would be read as three symbols.

## Step 3. The checks against the sources

In the Boyer-Lindquist chart the check holds, before anything is written:
the chart at $\ell = 0$ to the published Boyer-Lindquist chart of `kerr`, with $GM/c^2 = m$;
the Ricci scalar to Simpson and Visser's (3.1), $R = 2\ell^2me^{-\ell/r}/(\Sigma r^3)$;
the Kretschmann scalar to their (3.3), with $\Xi = \ell\Sigma/2r^3$;
the mixed Einstein tensor $G^\mu{}_\nu$ to having the characteristic polynomial $(\lambda - \lambda_1)^2(\lambda - \lambda_2)^2$ with $\lambda_1 = -2\ell me^{-\ell/r}/\Sigma^2$ and $\lambda_2 = -2\ell me^{-\ell/r}(\Xi - 1)/\Sigma^2$, the eigenvalues of their orthonormal (3.11), so the source is of Hawking-Ellis type I with $\rho = -p_r$;
and the chart at $a = 0$ to their static (2.9), $f = 1 - 2me^{-\ell/r}/r$.
With $G$ and $c$ restored, (3.1) and (3.3) keep their form, since $m$, $a$ and $\ell$ are already lengths.

## Step 4. How the values are written

`EyeOfTheStormForms` writes every value in $\Sigma$ and $\Delta$.
It holds each power $e^{-k\ell/r}$ as a power of one symbol, so that every value is a rational function, factors it as Kerr-de Sitter's values are factored, writes each factor that is a multiple of $\Sigma$ or $\Delta$ by its name, and writes the exponential back.
$g_{tt}$, $g_{t\phi}$, $g_{\phi\phi}$ and the inverse components of the Boyer-Lindquist chart are written by hand in Kerr's grouping and read back against the computed values.
