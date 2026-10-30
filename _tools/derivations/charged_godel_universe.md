# Charged universes of Gödel's kind

## Sources and scope

The founding record is A Banerjee and S Banerji, Journal of Physics A: General Physics 1 (1968), 188-193, DOI 10.1088/0305-4470/1/2/302.
Its full text has not been read.
The readable equation source is P Klepáč and J Horský, arXiv:gr-qc/0009022, section 3.2.4, equations (3.9)-(3.10) at C=0.
That section explicitly identifies the homogeneous solution with equations (11)-(12) of the founding work.
The journal title is Charged universes of the Gödel type with closed timelike curves, Classical and Quantum Gravity 17 (2000), 2547-2555, DOI 10.1088/0264-9381/17/13/308.
Firstmate approved this equation source on 2026-10-30 after the founding-paper access routes were exhausted.
No quotation from the founding work is used.

## Homogeneous metric and source split

Use signature (-+++) and Gaussian geometrized units G=c=1 in this derivation.
Write ds²=-(dt+H dφ)²+dr²+D²dφ²+dz², with H'/D=2w and D''/D=m².
The positive-curvature-parameter branch has D=sinh(mr)/m and H=2w(cosh(mr)-1)/m².
In the orthonormal coframe (dt+H dφ,dr,D dφ,dz), the mixed Einstein tensor is diag(m²-3w²,w²,w²,m²-w²).
A magnetic field F=B D dr∧dφ has 8π times its mixed stress tensor equal to B² diag(-1,1,1,-1).
Taking Λ=-m²/2 and B²=w²-m²/2 leaves dust with 8πρ=2w².
Maxwell's equation ∇νF^{μν}=4πJ^μ gives 4πJ^μ=(2wB,0,0,0), and F_{μν}J^ν=0.
These statements have been checked independently with the collection's Geometry on an equivalent local exponential chart.
The eventual chart checks must verify them numerically on every displayed chart.

## Circular timelike curves

A circle at constant t,r,z has squared tangent norm D²-H².
For the hyperbolic branch this is (4/m²)sinh²(mr/2)[1+(1-4w²/m²)sinh²(mr/2)].
It becomes negative beyond r_c=(2/m)atanh(m/(2w)), whenever 0<m<2w.
The pressureless magnetic solution requires m²≤2w² and therefore has circular timelike curves.
The uncharged endpoint B=0 is Gödel's dust universe.
The m=0 member is the homogeneous Som-Raychaudhuri solution already in the collection.
The present page must show the dependence on m and B explicitly so it describes the larger charged family.

## Remaining work

History approval, complete chart definitions, per-chart numerical checks, stored curvature, every required drawing, related links, bibliography and page checks remain to be completed.
The trigonometric branch and the electromagnetic field attributed to Rebouças must be derived independently before any of their signs or numbers is published.

## Euclidean embedding and the axis

The plane of constant t and z has circular radius rho=sqrt(D²-H²) and embedding height derivative squared 1-rho'².
The Euclidean endpoint precedes the first null circle on every branch drawn.
At ω=m=1 its radius is 2 asinh(sqrt((sqrt(3)-1)/3))=0.9516430642, before the null circle log(3)=1.098612289.
At ω=1 in the linear member the endpoint is sqrt(3)/2, before the null circle 1.
At ω=k=1 in the trigonometric member the endpoint is 2 asin(sqrt((3-sqrt(5))/5))=0.8031727833, before the null circle 2 atan(1/2)=0.9272952180.
The trigonometric metric component printed as 5 sin²(r)+8 cos(r)-8 loses digits near its regular axis.
With s=sin²(r/2), exactly the same component is 4s(1-5s), and the embedding defect is 4s(4-30s+25s²)/(1-5s).
The defect tends to 4r² at the axis; its apparent 0/0 in the first expression is numerical cancellation, not a geometric singularity.
The existing half_angles helper rewrites both expressions before quadrature, so the integrand is evaluated continuously at the axis without sampling a singular expression.
No shared Slice implementation is changed and neither endpoint is entered as a horizon in Slice.exact.
The local finite-endpoint slope uses rho'=1 at the axis and a real height derivative tending to zero at the numerically located outer endpoint.

## Drawn domains and curvature

The spacetime diagrams unroll the cylinder of t and φ at half and three halves of the first null radius on each branch.
Their null slopes are c dt/dφ=-H-D and D-H; the curve moving right at half that radius is a circular null geodesic.
The diagrams identify their angular edges and distinguish a closed circular timelike curve from a freely falling light ray.
Every drawing fixes ω=1; the hyperbolic branch also fixes m=1 and the trigonometric branch k=1.
The trigonometric chart covers the first open polar interval 0<r<π/k; the axis is regular as a limiting point of its spatial embedding.
The curvature generator reads both scalars from the stored chart, with no curvature recomputation.
It draws each scalar against proper radius on 0≤ωr≤2, separately on each branch, and checks spatial constancy at radii on both sides of the first null circle.
The resulting Kretschmann values are 24, 44 and 72 in units of ω⁴, and the Ricci values are 0, 2 and 4 in units of ω², for the hyperbolic, linear and trigonometric drawing parameters respectively.
These constant scalars show that the first null circle is not a curvature singularity.
The conformal generator records this spacetime in NOT_DRAWN because a global time-oriented radial conformal square does not encode its closed angular timelike curves.
All views are stationary snapshots, so no movie or stack registry entry is required.
The shared view_text count is below 400 for every view.
The inside cylinder on each branch marks its own embedding at t=0.
The outside cylinders lie beyond the embedded spatial region and carry no such mark.
The slice exceptions record both that geometric exclusion and the different curvature parameters of the other branches, so moments from different universes are not overlaid.
