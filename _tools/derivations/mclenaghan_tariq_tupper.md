# McLenaghan-Tariq-Tupper electrovacuum

## Sources and access

Raymond G. McLenaghan and Nessim Tariq published their solution in Journal of Mathematical Physics 16, 2306-2312 (November 1975), DOI 10.1063/1.522461.
N. Tariq and B. O. J. Tupper published a class of algebraically general non-null electrovacua in General Relativity and Gravitation 6, 345-360 (August 1975), DOI 10.1007/BF00761969.
Both original full texts remain unread; their primary bibliographic records were fetched.
Tupper's sequel, General Relativity and Gravitation 7, 479-486 (June 1976), DOI 10.1007/BF00766405, acknowledges the independent McLenaghan-Tariq solution on its publisher record.
No quotation or equation is attributed to an unread original.

Alexey Anatolievich Magazev's 2021 article, DOI 10.1007/s11040-021-09385-3, is available as arXiv:1908.10035v3.
Section 5, equation (29), supplies an explicit chart, equations (31)-(33) its invariant frame and coframe, and equation (42) the integrated Hamilton-Jacobi equation.
The entire preprint was read.
Anatoly Leonidovich Evtyagin and Alexey Anatolievich Magazev's 2021 article, DOI 10.24147/1812-3996.2021.26(1).21-26, was read in the university-hosted journal issue, pages 21-26.
Its Section 2 states the difference between the metric and electromagnetic symmetries, and Sections 3-4 treat the Klein-Gordon equation.

Ian Anderson and Charles Torre's 2022 article, DOI 10.1007/s10714-022-02913-8, was read as arXiv:2109.02558v3.
Its introduction identifies the McLenaghan-Tariq electrovacuum as homogeneous and non-inheriting, and equations (2.4)-(2.14) explain duality rotations and symmetry inheritance.
Its new spacetime has a positive cosmological constant and is a distinct solution.
Its global hyperbolicity and geodesic completeness results are not transferred to this electrovacuum.

## Chart definitions

Magazev's equation (29), with the overall sign changed to the site's convention, is

\[
ds^2=-(dX^1+dX^2-kX^3dX^4)^2+(dX^2)^2+e^{-kX^2}(dX^3)^2+e^{kX^2}(dX^4)^2,
\qquad k>0.
\]

The substitution

\[
ct=X^1+X^2,\quad r=X^2,\quad Y=X^3,\quad Z=X^4,\quad a=2/k
\]

produces the exponential form

\[
ds^2=-(c\,dt-2Y\,dZ/a)^2+dr^2+e^{-2r/a}dY^2+e^{2r/a}dZ^2.
\]

Then

\[
x=a e^{r/a},\quad y=Y,\quad \varphi=Z/a
\]

produces the upper-half-plane form

\[
ds^2=-(c\,dt-2y\,d\varphi)^2+\frac{a^2}{x^2}(dx^2+dy^2)+x^2d\varphi^2.
\]

This last form is also explicitly displayed by Frausto, Kolář, Málek and Torre (2025), equation (4.16), with attribution to the textbook's equation (12.21).
The symbol \(\varphi\) is a real coordinate; no periodic identification is imposed by these transformations.
These changes of coordinates have yet to be checked by the page's generator.

## Pending verification

Determine the Maxwell two-form and its inherited subgroup independently of the metric Killing fields.
Compute and print the chart tensors with the existing generator.
Check the apparent coordinate inconsistency in the extracted text of Evtyagin and Magazev's equation (1) before using that matrix.
Record every diagram's slice and any physically justified absence before drawing.
History publication awaits the keyed review decision.
