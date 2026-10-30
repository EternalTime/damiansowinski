# Kerr's black hole in a universe, and Kerr radiating

## Sources and scope

P C Vaidya and L K Patel, Radiating Kerr Metric, Physical Review D 7 (1973), 3590-3593, DOI 10.1103/PhysRevD.7.3590.
Its body sits behind the publisher's paywall and has not been read; its abstract has.
The radiating metric is taken from the authors' own restatement, P C Vaidya and L K Patel, Radiating Kerr particle in Einstein universe, Pramana 32 (1989), 731-739, section 2, equations (1)-(5), read in full.
P C Vaidya, The Kerr metric in cosmological background, Pramana 8 (1977), 512-517, is read in full from the Indian Academy of Sciences scan; its equations (13)-(15), its p. 515 and its section 4 are the source of the other two spacetimes.
Several citing papers print the 1977 pages as 151; those are the pages of Vaidya's 1984 Pramana paper on the de Sitter background.

The page holds three spacetimes, listed in `REGIONS` in `metric_tags.py`: Kerr radiating, Kerr in Einstein's static universe, and Kerr in an expanding universe.
Units are c = G = 1 and the signature is (-,+,+,+); Vaidya and Patel write (+,-,-,-), so every line element is theirs with its sign changed, which leaves the Ricci tensor as it is.

## Kerr radiating

In the barred chart of the 1989 paper, (4) and (5), with r = -q = T - u - (a/2) sin^2 alpha,
ds^2 = -2 theta (dT + (a q/2A) sin^2 alpha d beta) + (A^2 cos^2 alpha + q^2)(d alpha^2 + sin^2 alpha d beta^2) + (1 + a^2 sin^2 alpha/4A^2 - 2mq/(A^2 cos^2 alpha + q^2)) theta^2,
with theta = du + A sin^2 alpha d beta, W = A sin(alpha), A^2 = k^2 + a u, and m = C/A^3 from the paper's M W^3 cosec^3 alpha = constant.
The paper writes N = V W_u + W_alpha in its first chart (1)-(2); read in the barred chart, the quantity that makes the metric a solution is N = dW/d(alpha) at fixed u-bar, A cos(alpha), and N^2 + q^2 = r^2 + A^2 cos^2 alpha.
The scan also prints the coefficient of theta^2 as 1 - W_u^2 - 2Mq/(N^2 + q^2); the field equations need 1 + W_u^2 - 2Mq/(N^2 + q^2), with W_u^2 = a^2 sin^2 alpha/4A^2.
With either literal reading the Ricci tensor is not that of pure radiation (R_uT = 4a^2 sin^2 alpha cos^2 alpha/Sigma^2 with the minus sign, at every mass), so the scan carries sign or chart slips in N and in W_u^2; the page uses the form that solves the field equations, and `vaidya_patel_kerr_check` in `print_charts.py` holds it.
At a = 0 both readings are Kerr's metric, so the slips touch only the radiation.

The checks, at thirty digits at three points of each radiating chart, with a = 3/10, k = 7/10 and C = 1/2:
the Ricci tensor is sigma theta_mu theta_nu with theta null, which is pure radiation along one null direction;
sigma is positive at every point checked, all of them at r > 0 with the mass falling;
with C = 0 every component of the Riemann tensor vanishes, so the spacetime is flat however A runs;
and the retarded chart is the barred chart pulled back along T = r + u + (a/2) sin^2 alpha.
At a = 0, A = k and m = C/k^3 are constant and the metric is Kerr's in Kerr-Schild form; `_tools/test_vaidya_patel_kerr.py` holds its Kretschmann scalar to Kerr's there.
For a > 0 the mass C/A^3 and the angular momentum mA = C/A^2 both fall as A grows.

A and m are held as functions of u with the slopes `vm.RATES` declares, A' = a/(2A) and m' = -3am/(2A^2), so every printed value is rational in them.

## Kerr in Einstein's static universe

Vaidya's (13) and (14): ds^2 = -2 theta dt + (1 + 2 m mu) theta^2 + Sigma (R^2 d alpha^2/(R^2 - a^2 sin^2 alpha) + sin^2 alpha d beta^2), theta = du + a sin^2 alpha d beta,
Sigma = (R^2 - a^2) sin^2(r/R) + a^2 cos^2 alpha, mu = R sin(r/R) cos(r/R)/Sigma, r = t - u.
His p. 515 gives a perfect fluid with 8 pi p = -(1 - 2 m mu)/R^2 + Lambda and 8 pi rho = 3(1 - 2 m mu)/R^2 - Lambda.
In this signature the mixed Einstein tensor then has the eigenvalue -(1 - 2 m mu)/R^2 three times and -3(1 - 2 m mu)/R^2 once, which the check holds at three points of each of the two charts, with m = 1/2, a = 2/5 and R = 3.
In Vaidya's chart u and t enter only through the radius r = t - u, which is held as a function with slopes -1 and 1, as A and m are in the radiating charts; with the sines of (t - u)/R written out sympy did not finish building the chart in ten minutes.

The Boyer-Lindquist chart is derived here.
Write Vaidya's chart in r = t - u, the form his section 4 locates the horizon in, put dt -> dt + f(r) dr and d beta -> d phi + h(r) dr, and ask that every cross term with dr vanish for every alpha.
The two conditions give h = a/Delta and f = -2 m R sin(r/R) cos(r/R)/Delta, with
Delta = (R^2 - a^2) sin^2(r/R) + a^2 - 2 m R sin(r/R) cos(r/R),
which depends on r alone because mu Sigma does; then g_rr = Sigma/Delta.
Delta = 0 divided by cos^2(r/R) is R^2 tan^2(r/R) - 2 m R tan(r/R) + a^2 = 0, Vaidya's (18).
As R -> infinity Delta is Kerr's r^2 - 2mr + a^2 and the chart is Boyer and Lindquist's with the sign of a reversed.
The check holds the chart to Vaidya's chart pulled back along that Jacobian, with u = t - r: his chart's components depend on u and t through r alone, so the event (t - r, alpha, phi, t) serves as the image of (t, r, alpha, phi) at every point checked.

## Kerr in an expanding universe

Vaidya's (15) is e^{2F(t)} times (13) with m e^{-2F} in place of m, and F is left free by the field equations.
The check declares F = 3t/10 + t^2/10, any F serving, and holds the mixed Einstein tensor at three points to real eigenvalues of which exactly one pair is equal: the two pressures across the direction theta, with a third along it, as Vaidya's p. 515 says.
## Charts left out

Vaidya's charts in r = t - u in place of u, for the static and the expanding universe, were printed and checked and then left out: each wrote about 0.7 MB of trigonometric polynomials, with Sigma and Delta nowhere factored by name, and the page would have passed four megabytes.
The Boyer-Lindquist chart carries the static universe in t and r, and Vaidya's own charts carry both universes in his coordinates.

## Conformal diagram

None is drawn, and `NOT_DRAWN` in `conformal.py` lists the spacetime.
The radiating spacetime is time dependent with no horizon known in closed form, so no chart of it is conformally compact in closed form.
The expanding universe leaves F free, so any diagram would draw one choice of F as if it were the spacetime.
The static universe's axis runs from the black hole through its two horizons to the antipode, and its equator ends on the naked antipodal ring, which no two dimensional diagram of the axis shows.

## Drawings

Every chart is drawn on its axis, alpha = 0, where the rotation about the axis keeps light rays launched in the plane of the time and the radius in that plane.
The universe is drawn at m = 1, a = 3/5 and R = 4, where R tan(r/R) = m +- sqrt(m^2 - a^2) puts the horizons at r = 1.6914 and 0.1998, inside Kerr's 1.8 and 0.2, out to r = pi R/2.
The expanding universe declares F = t/10, under which Vaidya's condition a^2/m^2 < e^{-4F} holds while t < (5/2) ln(25/9) = 2.554.
The radiating hole is drawn at k = 3/5, C = 27/125 and a = 1/10, so its mass is 1 at u = 0; its views stop at r = 0, since beyond the disc u + r stops being a time function on the axis.
The embedding is the equator of the moment t = 0 of the Boyer-Lindquist chart, through the throat on the outer horizon, of circumference radius 2mR cos(r_+/R)/sqrt(R^2 - a^2) = 1.8447, to the level circle r = 9.741, beyond which no surface of revolution in flat space carries it; nearer the antipode the circles turn timelike.
That moment is marked on the Boyer-Lindquist axis; `exceptions/vaidya_patel_kerr.json` says why the other views do not mark it.
The curvature figure maps K and R of the Boyer-Lindquist chart in the oblate frame about the ring, out to 4m, with R = 6(1 - 2m mu)/R^2 checked against the stored scalar; R turns negative near the ring, where Vaidya's density does.
No view steps through time, so no entry in the tables of movies or stacks is needed.
