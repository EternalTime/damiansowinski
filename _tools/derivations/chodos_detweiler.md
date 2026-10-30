# A five-dimensional universe with compact spatial circles

## Sources and charts

Chodos and Detweiler, Physical Review D 21, 2167-2170 (1980), DOI 10.1103/PhysRevD.21.2167, specify four periodic spatial coordinates with the same coordinate period $L$.
Their unnumbered metric at the top left of page 2168 specializes (2) and (3) to exponents $(1/2,1/2,1/2,-1/2)$.
Chyba, American Journal of Physics 53, 863-872 (1985), DOI 10.1119/1.14353, restates the metric in (6.5) and its rescaled form in (6.6), using the opposite signature.

The founding article's publication metadata was verified against the publisher record and Crossref.
Its readable cached extraction was consulted in the research conversation but was not archived as a full-text file, and subsequent retrieval is currently unreproducible.
A later saved HTTP 502 response contains no original-paper text.
The retained institutional author scan of Chyba independently supports the metric, identifications and circumference normalization in (6.2)-(6.7), page 868.

The two proposed charts retain all five dimensions:

- `synchronous`: $ds^2=-c^2dt^2+(t/t_0)(dx^2+dy^2+dz^2)+(t_0/t)dw^2$, with every spatial coordinate periodic with period $L$.
- `rescaled`: $ds^2=-c^2dt^2+(t/\tau)(dX^2+dY^2+dZ^2)+(\tau/t)dW^2$, with the coordinate periods transformed with the coordinates.

Here $t>0$, $t_0>0$ and $\tau>0$ have units of time, and $L$ and all spatial coordinates have units of length.
The source's indices $1,2,3,5$ are written $x,y,z,w$.
The positive-time branch is the cosmology described by the authors.

## Coordinate normalization

For a fixed reference time $\tau$, define

$$X=\sqrt{\tau/t_0}\,x,\quad Y=\sqrt{\tau/t_0}\,y,\quad Z=\sqrt{\tau/t_0}\,z,\quad W=\sqrt{t_0/\tau}\,w.$$

The time coordinate is unchanged.
Substitution in the original line element produces $t/\tau$ on the first three spatial terms and $\tau/t$ on the fifth-direction term.
The periods become $L\sqrt{\tau/t_0}$ for $X,Y,Z$ and $L\sqrt{t_0/\tau}$ for $W$.
The coefficient $\tau/t$ agrees with Chyba's (6.6), inspected on page 868 of the author-hosted scan.
The fifth period agrees with his (6.7), with $2\pi r_5=L$.
The cached extraction of the 1980 paper transcribes this coefficient incorrectly as $t/\tau$.

The proper circumferences on a surface of cosmic time are

$$C_x=C_y=C_z=L\sqrt{t/t_0},\qquad C_w=L\sqrt{t_0/t}.$$

Both charts describe these same lengths.
At $t=\tau$, the rescaled spatial metric is the identity, while its four coordinate periods differ.
At $t=t_0$, the four circles have the same proper circumference $L$, while their expansion rates have opposite signs.
Neither statement requires a change in the number of spacetime dimensions.

Each ratio in a metric coefficient is dimensionless; each coordinate period has units of length.
The chart for tensor computation uses $T=ct$, so $g_{TT}=-1$ and the spatial coefficients are $T/(ct_0)$ and $ct_0/T$.

## Vacuum equations and scale factors

The scale factors are $a(t)=\sqrt{t/t_0}$ and $b(t)=\sqrt{t_0/t}$.
Their logarithmic rates are $\dot a/a=1/(2t)$ and $\dot b/b=-1/(2t)$.
The four Kasner exponents sum to one and their squares sum to one.
The spatial volume of one fundamental cell is $L^4a^3b=L^4t/t_0$.

The source normalization and its exact coordinate pullback were checked in `/Users/owlshome/BootLoops/bootloops/.venv/bin/python` using the site's `verify_metrics.Geometry`.
All 25 components of the Ricci tensor vanish for each five-dimensional chart.
The two metric matrices and all four proper circumferences agree at three rational points with distinct values of $t_0$ and $\tau$.
The task's `normalization-check.py` and its log preserve the calculation; its timed computation took 0.063 seconds.
The generator's chart checks will carry these source identities when the page is implemented after review.

## Drawings

Drawing generation follows the chart implementation and its source checks.
