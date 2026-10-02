# Mass inflation

The three charts of `mass_inflation.json` are written by `print_charts.py --metric mass_inflation`, which computes every tensor from each line element with the checker's own `Geometry` and reads every printed value back before it writes.
This note records the source of each chart, the field equations each is held to, and what the diagrams draw.

## Step 1. The spacetime

Eric Poisson and Werner Israel take a charged, spherical black hole crossed by two radial streams of lightlike particles that do not interact, one falling in and one running outward (`poisson1989`).
Their metric is $ds^2 = g_{ab}dx^adx^b + r^2d\Omega^2$, and their (1) defines the mass function,

$$1 - \frac{2m}{r} + \frac{e^2}{r^2} = f = g^{ab}(\partial_a r)(\partial_b r),$$

with $e$ the charge, constant.
The page writes the charge as the charge radius $r_q$ of the Reissner-Nordström page, so that the letter $q$ stays free for the null coordinate of the conformal diagram, and takes $G = c = 1$, as every source does.

## Step 2. The ingoing chart

With the influx alone the metric is the charged Vaidya metric, which Poisson and Israel write above their (6) as $ds^2 = 2\,dr\,dv - f\,dv^2 + r^2d\Omega^2$ with $m = m_L(v)$, citing Sullivan and Israel, and which is Bonnor and Vaidya's metric with the charge held constant (`bonnor1970`).
Hiscock used it in 1981 for the inside of the charged hole (`hiscock1981charged`, as `brady1999` section 3.2 reports it).
`mass_inflation_check` holds the chart to $G^a{}_b = \mathrm{diag}(-1, -1, 1, 1)\,r_q^2/r^4$ plus the one component $G_{vv} = 2m'(v)/r^2$ of the null dust, and, at constant mass, to the published metric of `rn_metric` with $r_s = 2m$, pulled back along $t = v - r_*$.

## Step 3. The double null chart

Where both streams flow Poisson and Israel write the plane of radial motion in null form, $-2e^{2\sigma}dU\,dV$, their text above (7), which Bonanno, Droz, Israel and Morsink repeat as their (20) (`bonanno1995`).
The chart leaves $r(u, v)$ and $\sigma(u, v)$ free, and no component assumes a field equation.
With $E_{ab}$ the Maxwell field of the charge, whose one component on the plane is $E_{uv} = r_q^2e^{2\sigma}/r^4$, the check holds, for any $r$ and $\sigma$,

$$\partial_v m = -\tfrac{1}{2}r^2e^{-2\sigma}\left(\partial_u r\,G_{vv} - \partial_v r\,(G_{uv} - E_{uv})\right),$$

and the same with $u$ and $v$ exchanged.
With $G_{vv} = 2L_{in}/r^2$ and $G_{uv} = E_{uv}$ this is $\partial_v m = -e^{-2\sigma}L_{in}\,\partial_u r$, the first of their (2).
Where the two equations with no flux in them hold, $G_{uv} = E_{uv}$ and $G_{\theta\theta} = E_{\theta\theta}$, it also holds

$$\partial_u\partial_v m = \tfrac{1}{4}r^3e^{-2\sigma}G_{uu}G_{vv} = \frac{e^{-2\sigma}L_{in}L_{out}}{r},$$

which is the second of their (3), $\Box m = -r^{-1}(4\pi r^2)^2T^{ab}T_{ab}$ with $\Box\psi = -2e^{-2\sigma}\psi_{,UV}$, and Bonanno's (44).
Its integral is their (7).
One stream alone gives a vanishing source, which is why the influx by itself leaves the mass bounded.

## Step 4. The advanced chart

Patrick Brady and John Smith write the general spherical line element as their (1), $ds^2 = -g\bar{g}\,dv^2 - 2g\,dv\,dr + r^2d\Omega^2$ (`brady1995`).
The reader takes no barred name, so the page writes $h$ for their $\bar{g}$ and says so where $h$ is described.
The check holds the chart, for any $g$ and $h$, to $G_{rr} = 2\,\partial_r g/(rg)$ and $G^v{}_v + G^r{}_r = 2(\partial_r(rh) - g)/(r^2g)$, which with a massless scalar field and the charge's Maxwell field are their (5) and (6), and to being the ingoing chart where $g = -1$ and $h = -f$.
Their (11) is the mass function, $m = \tfrac{1}{2}r(1 + e^2/r^2 - \bar{g}/g)$.

## Step 5. Ori's shell

Amos Ori keeps the influx continuous and makes the outflux one thin shell, with a charged Vaidya metric on each side (`ori1991inflation`); the equations here are as Bonanno, Droz, Israel and Morsink restate them, their (3) to (9), since Ori's letter was not to hand.
`ori_shell.py` integrates the model, and its opening lines derive the form it uses: along the shell, with $\lambda$ an affine parameter and $z_i = R/(dv_i/d\lambda)$, each side's geodesic equation gives $dz_i/d\lambda = (1 - r_q^2/R^2)/2$, so $z_2 - z_1$ is a constant $Z$, and $R' = f_iR/2z_i$ gives $m_2 = m_1 - ZR'$.
Raychaudhuri's equation on each side then agrees, so the influx is continuous across the shell, and Bonanno's (7), $dm_2/f_2 = dm_1/f_1$, follows; `consistent()` measures it on the numbers.
The drawings take $r_q = 0.96\,m_0$, the charge of Reissner-Nordström's own diagrams, the tail $m_1 = m_0 - (m_0/50)(v_0/v)^{11}$ from $v_0 = 10\,m_0$ on and $0.98\,m_0$ before, Price's $p = 12$, and a shell of mass $m_0/50$ that crosses $v_0$ at $r = m_0$.
Behind the shell the advanced time reaches the Cauchy horizon at a finite value, taken as zero, and $m_2$ grows as $w^{-12}e^{\kappa w}$ in the advanced time $w$ before the shell, Bonanno's (9): it passes $2\,m_0$ at $5 \times 10^{-11}\,m_0$ from the Cauchy horizon.
An outgoing ray behind the shell still arrives there at a radius above zero, since $m_2$ is integrable in its own advanced time, which is the weakness Ori found.
Every function of the shell has a jump in its derivative on the ingoing ray $v_0$, where the influx is switched on, so each spline and each integration is cut there.

## Step 6. The diagrams

The spacetime diagrams draw the ingoing chart twice: the tail falling in, with the event horizon and the outgoing ray the shell runs along marked, and the chart behind the shell, hatched beyond the shell and beyond the Cauchy horizon $v = 0$.
`--verify` checks the outgoing rays of both against `ori_shell.label_before` and `label_behind`, integrators the tracing does not use.
The conformal diagram draws Ori's model from the ingoing ray $10\,m_0$ before $v_0$: an ingoing ray is named on both sides by its advanced time $w$ before the shell and an outgoing ray by where it began, as the docstring of `mass_inflation` in `conformal.py` states.
The embedding diagram is the tail falling in alone, slices of constant $v - r$, a movie: the mass grows by a fiftieth, so it is Reissner-Nordström's funnel sinking by a little.
The inflation itself is in none of the three to the eye, since behind the shell it happens within $10^{-10}\,m_0$ of the Cauchy horizon; the caption of the view behind the shell gives the numbers.
The double null and advanced charts leave two functions free, which only a numerical evolution fixes, so they have the conformal diagram and no spacetime diagram of their own.

## Step 7. Names

Every first name printed is its paper's own author line or record: Michael Simpson and Roger Penrose (Crossref, `10.1007/BF00792069`), Patrick R. Brady and John D. Smith (the author line of arXiv gr-qc/9506067), William A. Hiscock, Richard H. Price, Tevian Dray and Gerard 't Hooft, Shahar Hod and Tsvi Piran, Mihalis Dafermos and Amos Ori (Crossref).
I. H. Redmount and J. S. F. Chan and K. C. K. Chan are printed with the initials their papers print; Robert Mann's first name is INSPIRE's.
Poisson and Israel's letter prints initials, and their names in full are the author line of `poisson1990`.
