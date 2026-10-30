# Relativistic self similar discs

## Sources

The founding journal scans are Lynden-Bell and Pineault, Monthly Notices of the Royal Astronomical Society 185 (1978), 679-694 and 695-712, DOI 10.1093/mnras/185.4.679 and 10.1093/mnras/185.4.695.
The original metric equations use c=1.

## Chart inventory

- Static axisymmetric Weyl coordinates (t,R,z,phi), Paper I equation (3.1), with vacuum potentials obeying equations (3.2)-(3.4) off the disc.
- Infinite counter rotating disc in spherical coordinates (t,r,theta,phi), Paper I equations (3.14)-(3.16).
- Counter-rotating massless-particle specialization, Paper I equation (4.2), at n=1/2.
- Finite counter rotating disc in Weyl or spherical coordinates, Paper I section 3.2, with potential series and a quadrature for the second potential.
- Rotating self similar coordinates (t,r,theta,phi), Paper II equation (2.2), specialized to k+n=1 by its thin pressureless-disc boundary conditions.

The finite family and the infinite family have different outer boundary conditions.
A finite-disc outer radius and the arbitrary scale used to normalize the infinite metric have different meanings.
The rotating angular functions require the boundary-value integration of Paper II sections 3.1 and 3.2.
Its WKB expressions in the Appendix are approximations and must not be used as an exact rotating metric.

## Matter and vacuum

The source is an infinitesimally thin equatorial matter surface.
Paper I writes a delta source at z=0 in equation (3.7), and Paper II integrates its field equations across theta=pi/2 on pages 701-702.
Smooth exterior tensors apply on either open side of this surface.
Vacuum Ricci components away from the disc do not describe its surface stress tensor.
For the infinite static family, the smooth curvature diverges at the origin for 0<n<1/2.
The massless member n=1/2 has a flat vacuum exterior and a nonzero distributional sheet.
The change rho=2 sqrt(Ar) sin(theta/2), zeta=2 sqrt(Ar) cos(theta/2) makes its spatial metric Euclidean and its lapse zeta/(2A).

## Parameters and limits

For the static infinite family, n=V^2/(1+V^2), where V is the speed of either stream measured by a fixed observer in units c=1.
The interval 0<n<1/2 describes massive streams, and n=1/2 describes counter rotating massless particles.
For the rotating family, V is the speed measured by a locally nonrotating observer.
Paper II obtains the end of its pressureless rotating sequence numerically at V approximately 0.438c.
The drawings use A=1, n=1/5 for the infinite counter rotating member and V=0.2c for the cold rotating member.
The rotating boundary solve gives n=0.045535001497421096 and Q'(0)=0.020898796121853757.
Its functions are integrated from the regular-axis expansions, with N(0)=Z(0)=Q(0)=0 and e^P=sin(theta).
The exact vacuum equations and cold-disc jump conditions determine the angular functions; no WKB expressions enter the calculation.

## Drawings

The constant-time, constant-polar-angle surfaces have conical Euclidean embeddings.
Paper I section 3.3 and Paper II page 707 state the constructions.
The source surface and vacuum sides must be identified in any drawing crossing the equator.
The embedded surfaces are the vacuum cones theta=pi/3 at t=0, over 0.3A<=r<=3.5A.
Each is a single moment, and each view has one surface.
Radial light views at theta=pi/3 and phi=0 show null curves of the induced metric, without a claim that the curves are spacetime geodesics.
The Weyl view shows null geodesics along the upper axis of the infinite counter rotating member.
The conformal view is restricted to the totally geodesic upper axial section of the static n=1/5 member.
Its r=0 boundary is timelike and singular; its diagram asserts no global causal classification for the complete axisymmetric spacetime or the rotating family.
The curvature view reads the stored static Kretschmann scalar along theta=pi/3; it contains no distributional surface scalar.

## Source comparison and verification

The retained source images were inspected at Paper I pages 683, 686, 687, 688 and 690, Paper II pages 696, 698 and 699, and Cai-Shu pages 4 and 9.
The finite-disc potential uses the even Legendre polynomials in the outer series on Paper I page 687; the earlier Newtonian expression on page 683 has a different printed subscript.
The finite family is represented by the Weyl chart with its harmonic potential, disc Dirichlet value, zero exterior normal derivative and asymptotic normalization.
The spherical chart is the explicit infinite member of that construction; a second finite spherical chart would repeat the coordinate map with the generic Weyl potentials.
In the cold rotating member k+n=1, Cai-Shu's angular variable Theta equals the original polar angle and their radial exponent is n/(1-n).
Their printed relation Z=Z_LP-ln(k) in equation (2-5) omits a factor of two: pulling back the radial line element gives Z=Z_LP-2ln(k).
The implementation uses the original radial chart and only angular derivatives of Z, so this additive constant does not enter it.
Each chart has pointwise source comparisons with no symbolic determinant or simplify call in its check.
The dedicated test checks the surface stress through the equatorial derivative jump and independently pulls the massless member back from Minkowski space.
Reconcile shared registrations and reciprocal relations only on the serialized landing turn.

## History

The gravity of an infinite disc can support circular motion at a speed independent of radius.
Donald Lynden-Bell and Serge Pineault constructed such relativistic equilibria in 1978, publishing their counter rotating and rotating families in Monthly Notices of the Royal Astronomical Society [lyndenbell1978, lyndenbellpineault1978b].
The constants $G$ and $c$, together with the circular speed, supply no characteristic length.
These infinite discs are self similar, with a change of radial scale accompanied by a change of the time scale [lyndenbellpineault1978b].

Equal streams of stars circle in opposite senses in the counter rotating family.
Each stream is cold, and its stars follow circular orbits with a speed independent of radius.
The opposing streams have equal and opposite angular momenta, and their gravitational field is static, with no dragging of inertial frames [lyndenbell1978].
Lynden-Bell and Pineault also constructed finite discs whose geometry approaches the infinite disc solution well inside the outer edge [lyndenbell1978].

The matter occupies an infinitesimally thin equatorial sheet.
Its surface energy density and azimuthal stress enter Einstein's equations through a delta function supported on that sheet [lyndenbell1978].
The regions above and below the disc satisfy the vacuum Einstein equations, with a vanishing smooth Ricci tensor away from the matter surface [lyndenbell1978, lyndenbellpineault1978b].

At the speed of light, the counter rotating streams consist of particles with zero rest mass [lyndenbell1978].
Lynden-Bell and Pineault obtained analytic photon trajectories for this limiting disc.
Some photons circle at constant radius while oscillating across the equatorial plane, and others spiral inward or outward [lyndenbell1978].

Matter rotating in one sense drags the surrounding inertial frames.
Lynden-Bell and Pineault reduced the stationary field equations through self similarity and determined the rotating disc's angular metric functions numerically [lyndenbellpineault1978b].
An ergoregion appears around the equatorial plane at a circular speed of approximately $0.384c$ measured by locally nonrotating observers, where the stationary Killing vector becomes spacelike.
Its conical boundary extends toward the symmetry axis as the rotation speed increases.
Their cold rotating sequence ends near $0.438c$, with a diverging frame dragging parameter and an ergoregion extending throughout space away from the axis [lyndenbellpineault1978b].

The infinite discs have infinite total mass and a geometry that is not asymptotically flat [lyndenbell1978, lyndenbellpineault1978b].
The cold counter rotating equilibria are unstable to gravitational fragmentation into rings [lyndenbell1978].
Mike J. Cai and Frank H. Shu included isothermal pressure in the disc plane in 2002, retaining infinitesimal vertical thickness and a circular speed independent of radius [caishu2002].
