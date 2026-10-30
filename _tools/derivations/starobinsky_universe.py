"""The homogeneous semiclassical geometry of Starobinsky (1980), Eqs. (4), (7), (12).

All times in the numerical solution are H_0 tau and all lengths are H_0^-1.
The finite initial displacement uses a second-order expansion at the de Sitter saddle.
The evolution retains the full homogeneous equation, including the anomaly term.
"""
from functools import lru_cache

import numpy as np
import sympy as sp

import chart_printer as cp
import verify_metrics as vm

ID = "starobinsky_universe"
CHART_IDS = ("cosmic", "conformal", "closed_de_sitter", "flat_de_sitter", "open_de_sitter")
ANGLES = r"\left(d\theta^2 + \sin^2\theta\,d\phi^2\right)"


def chart_spec(system):
    """Two evolving FLRW charts and the three exact special solutions of Eq. (7)."""
    evolving = system in ("cosmic", "conformal")
    first = r"\eta" if system == "conformal" else r"\tau"
    radial = r"\chi" if system in ("closed_de_sitter", "open_de_sitter") else "r"
    coords = [first, radial, r"\theta", r"\phi"]
    if evolving:
        parameters = [f"a = a({first})", "k", "H_0", "M"]
        space = r"\dfrac{dr^2}{1-kr^2} + r^2" + ANGLES
        line = (r"ds^2 = a^2\left(-d\eta^2 + " + space + r"\right)" if system == "conformal"
                else r"ds^2 = -d\tau^2 + a^2\left(" + space + r"\right)")
        domains = [first + r"\in I", r"a>0", r"r\in[0,\infty)\;\text{for}\;k\le0",
                   r"r\in[0,1/\sqrt{k})\;\text{for}\;k>0"]
    else:
        parameters = ["H_0"]
        if system == "flat_de_sitter":
            space = r"dr^2 + r^2" + ANGLES
            factor = r"e^{2H_0\tau}"
        else:
            trig = r"\sin" if system == "closed_de_sitter" else r"\sinh"
            scale = r"\cosh" if system == "closed_de_sitter" else r"\sinh"
            space = r"d\chi^2 + " + trig + r"^2\chi" + ANGLES
            factor = scale + r"^2(H_0\tau)/H_0^2"
        line = r"ds^2 = -d\tau^2 + " + factor + r"\left(" + space + r"\right)"
        domains = [r"\tau\in(0,\infty)" if system == "open_de_sitter" else r"\tau\in(-\infty,\infty)",
                   (r"\chi\in[0,\pi]" if system == "closed_de_sitter" else radial + r"\in[0,\infty)")]
    domains += [r"\theta\in[0,\pi]", r"\phi\in[0,2\pi)"]
    domains = [domain.replace(r"\in", r" \in ") for domain in domains]
    probe = vm.Reader(coords, parameters, ())
    x = probe.symbol[first]
    pretty = cp.hyperbolic(probe.parameters["H_0"] * x) if system in ("closed_de_sitter", "open_de_sitter") else sp.factor
    return {"metric_id": ID,
            "system": {"id": system, "name": {"cosmic": "Cosmic Time", "conformal": "Conformal Time",
                       "closed_de_sitter": "Closed de Sitter", "flat_de_sitter": "Flat de Sitter",
                       "open_de_sitter": "Open de Sitter"}[system], "coords": coords, "domains": domains,
                       "parameters": parameters, "line_element": line},
            "chart_line_element": line, "printer": {"primed": ["a"]}, "pretty": pretty,
            "check": lambda chart: chart_check(chart, system)}


def chart_check(chart, system):
    """Check the independent FLRW Ricci formula at regular points, with no symbolic determinant."""
    x, r, theta, phi = chart.symbols
    p = chart.reader.parameters
    g = chart.geo.g
    ric = chart.geo.ricci_ll()
    evolving = system in ("cosmic", "conformal")
    for time, radius, curvature in ((sp.Rational(2, 5), sp.Rational(1, 3), -1),
                                    (sp.Rational(7, 5), sp.Rational(1, 4), 1)):
        point = {x: time, r: radius, theta: sp.Rational(4, 5), phi: sp.Rational(1, 5),
                 p["H_0"]: sp.Rational(3, 2), chart.reader.c: 1}
        if evolving:
            a = p["a"]
            values = {a: 2, sp.diff(a, x): sp.Rational(3, 7), sp.diff(a, x, 2): sp.Rational(-1, 9)}
            point.update({p["k"]: curvature, p["M"]: sp.Rational(1, 2)})
            A, B, C = (values[a], values[sp.diff(a, x)], values[sp.diff(a, x, 2)])
            if system == "cosmic":
                temporal, spatial = -3*C/A, C/A+2*(B*B+curvature)/A**2
            else:
                temporal, spatial = -3*(C/A-B*B/A**2), (C/A+B*B/A**2+2*curvature)/A**2
        else:
            values = {}
            temporal, spatial = -3*p["H_0"]**2, 3*p["H_0"]**2
        for i in range(4):
            for j in range(4):
                expected = temporal if i == j == 0 else spatial*g[i,j] if i == j else 0
                delta = (sp.sympify(vm._at(ric, (i,j))) - expected).subs(values).subs(point)
                if abs(complex(sp.N(delta, 30))) > 1e-24:
                    raise AssertionError(f"{ID}/{system}: FLRW Ricci component {i}{j} at {time}")


MASS = 0.5
DISPLACEMENT = 1e-4
EXPONENT = np.sqrt(2.25 + MASS**2) - 1.5


def trace_rhs(_, state, mass=MASS):
    """Regular trace equation, a'/a=h, h'=R/6-2h^2, R'=S, for k=0 and H_0=1."""
    N, h, R, S, eta = state
    Q = h*h
    return [h, R/6-2*Q, S, -3*h*S-mass**2*R+12*mass**2*Q*(R/6-Q), np.exp(-N)]


def constraint(h, R, S, mass=MASS, curvature=0.0):
    """M^2 times Eq. (4), using R=6(h'+2h^2+k/a^2), with H_0=1."""
    Q = h*h+curvature
    return mass**2*(Q-Q*Q)+Q*R/3-R*R/36+h*S/3


class Trajectory:
    """Full trace evolution with the exact constraint at the finite initial point.

    Time and a are normalized at h=1/sqrt(2): tau=0 and a=1 there.
    The second-order unstable-manifold expansion is the only initialization approximation.
    The sample extends 160 H_0^-1 beyond that crossing, through scalaron oscillations.
    """

    def __init__(self, displacement=DISPLACEMENT, tolerance=2e-12):
        from scipy.integrate import solve_ivp
        self.displacement = displacement
        q = EXPONENT
        b = q*(4*q+3)/(6*(q+1))
        h = 1-displacement
        hp = -q*displacement+b*displacement**2
        R = 6*(hp+2*h*h)
        S = (R*R/12-h*h*R-3*MASS**2*(h*h-h**4))/h
        initial = [0.0, h, R, S, 0.0]

        def crossing(_, state):
            return state[1]-1/np.sqrt(2)
        crossing.terminal = True
        first = solve_ivp(trace_rhs, (0, 250), initial, method="DOP853", events=crossing,
                          dense_output=True, rtol=tolerance, atol=tolerance*0.01, max_step=0.25)
        if not first.success or not len(first.t_events[0]):
            raise RuntimeError("Starobinsky trajectory did not reach h=1/sqrt(2)")
        self.shift = first.t_events[0][0]
        middle = first.y_events[0][0]
        self.N0, self.eta0 = middle[0], middle[4]
        self.scale = np.exp(self.N0)
        # Anchor eta at the crossing: subtracting two early-time conformal integrals
        # would lose all late-time digits after the long near-de Sitter expansion.
        self.past_eta = solve_ivp(lambda t, _: [np.exp(self.N0-first.sol(t)[0])],
                                  (self.shift, 0), [0.0], method="DOP853", dense_output=True,
                                  rtol=tolerance, atol=tolerance*0.01, max_step=0.2)
        normalized = [0.0, middle[1], middle[2], middle[3], 0.0]
        last = solve_ivp(trace_rhs, (0, 160), normalized, method="DOP853",
                         dense_output=True, rtol=tolerance, atol=tolerance*0.01, max_step=0.2)
        if not last.success:
            raise RuntimeError(last.message)
        self.first, self.last = first, last
        self.lo, self.hi = -self.shift, 160.0

    def state(self, time):
        time = np.asarray(time, dtype=float)
        if np.any((time < self.lo) | (time > self.hi)):
            raise ValueError(f"trajectory time outside [{self.lo}, {self.hi}]")
        shape = time.shape
        native = (time+self.shift).ravel()
        left = self.first.sol(np.minimum(native, self.shift))
        left[0] -= self.N0
        left[4] = self.past_eta.sol(np.minimum(native, self.shift))[0]
        right = self.last.sol(np.maximum(time.ravel(), 0))
        state = np.where(native <= self.shift, left, right)
        return state.reshape((5,)+shape)

    def a(self, time, order=0):
        N, h, R, _, _ = self.state(time)
        a = np.exp(N)
        return a if order == 0 else a*h if order == 1 else a*(R/6-h*h)

    def eta(self, time):
        """Conformal time from quadrature of 1/a in the same evolution."""
        return self.state(time)[4]

    def time(self, eta):
        from scipy.interpolate import CubicHermiteSpline
        if not hasattr(self, "inverse_eta"):
            times = np.linspace(-12, self.hi, 43001)
            etas = self.eta(times)
            self.inverse_eta = CubicHermiteSpline(etas, times, self.a(times), extrapolate=False)
        values = self.inverse_eta(np.asarray(eta, dtype=float))
        if np.any(~np.isfinite(values)):
            raise ValueError("conformal-time sample outside the numerical trajectory")
        return values

    def conformal_a(self, eta, order=0):
        time = self.time(eta)
        a, h, R = self.a(time), self.state(time)[1], self.state(time)[2]
        return a if order == 0 else a*a*h if order == 1 else a**3*R/6


@lru_cache(maxsize=1)
def trajectory():
    return Trajectory()


def declared_functions():
    """Scale factor and two derivatives for the existing numeric-function reader.

    evalf evaluates the stored scalar formula on the same double precision numerical
    trajectory; its digits are not additional precision of the differential equation.
    """
    made = {}
    for name, conformal in (("sta_a", False), ("sta_c", True)):
        for order in (2, 1, 0):
            label = name if order == 0 else name+str(order)

            def value(x, order=order, conformal=conformal):
                sol = trajectory()
                return sol.conformal_a(x, order) if conformal else sol.a(x, order)

            def evaluate(self, prec, value=value):
                if self.args[0].is_number:
                    return sp.Float(float(value(float(self.args[0]))), precision=prec)
            body = {"nargs": 1, "is_real": True, "_imp_": staticmethod(value), "_eval_evalf": evaluate}
            if order < 2:
                body["fdiff"] = lambda self, argindex=1, up=name+str(order+1): made[up](self.args[0])
            made[label] = type(label, (sp.Function,), body)
    return made


NUMERIC_INPUT = (r"Full semiclassical evolution; $M/H_0=1/2$, $a=1$ and $\tau=0$ at $h=H_0/\sqrt2$. "
                 r"Second-order saddle initialization at $h/H_0=1-10^{-4}$.")


def conformal_time(system, time):
    """Radial conformal coordinate in units H_0=1, with exact de Sitter primitives."""
    time = np.asarray(time, dtype=float)
    if system == "cosmic":
        return trajectory().eta(time)
    if system == "conformal":
        return time
    if system == "closed_de_sitter":
        return np.arctan(np.sinh(time))
    if system == "flat_de_sitter":
        return -np.exp(-time)
    return np.log(np.tanh(time/2))


def register_null_rays(nr):
    """Page-local rows for the existing ray tracer and its independent null integrals."""
    nr.DECLARED_FUNCTIONS.update(declared_functions())
    for system in CHART_IDS:
        radial = r"\chi" if system in ("closed_de_sitter", "open_de_sitter") else "r"
        first = r"\eta" if system == "conformal" else r"\tau"
        params = {"H_0": 1}
        functions = {}
        if system in ("cosmic", "conformal"):
            params.update({"k": 0, "M": MASS})
            functions = {"a": "sta_c(eta)" if system == "conformal" else "sta_a(tau)"}
        box = {"cosmic": (0,4,-2,2), "conformal": (0,3,-2,1),
               "closed_de_sitter": (0,np.pi,-1.5,1.5), "flat_de_sitter": (0,4,-2,2),
               "open_de_sitter": (0,3,0.2,3.2)}[system]
        nr.DIAGRAMS.append(nr.Diagram(ID, system, "radial", "$"+first+"$ and $"+radial+"$",
                           (first,radial),box,"$H_0r$" if radial == "r" else r"$\chi$",
                           "$H_0"+first+"$",params,nr.EQUATOR,tau=first.replace("\\", ""),
                           functions=functions,input=NUMERIC_INPUT if functions else None))
        nr.CAPTIONS[(ID,system,"radial")] = [
            "The radial plane ($\\theta=\\pi/2$, $\\phi=0$), with "
            + ("the evolving scale factor." if functions else "an exact de Sitter scale factor.")
            + " Light rays conserve $" + (r"H_0\eta" if radial != "r" else r"\eta") + r"\pm " + radial + "$."]
        nr.CLOSED_FORMS[(ID,system,"radial")] = (
            lambda t,r,s=system: conformal_time(s,t)+r,
            lambda t,r,s=system: conformal_time(s,t)-r,None)


def register_slices(sl):
    """The reference moment of each of the two embedded members of the family."""
    for system in ("cosmic","conformal","closed_de_sitter"):
        curved = system == "closed_de_sitter"
        view,source,coordinate = ("closed",system,r"\chi") if curved else ("flat","cosmic","r")
        sl.FLAT[(ID,system,"radial")] = (
            lambda view=view,source=source,coordinate=coordinate:
            sl.one(ID,lambda m:sl.along(0.0,*m.reach(source,coordinate)),view_id=view))
    sl.FLAT_METRICS.add(ID)
    for system in ("flat_de_sitter","open_de_sitter"):
        sl.HIDDEN[(ID,system,"radial")] = "a different member of the semiclassical family from the two embedded slices"


CONFORMAL_CAPTIONS = {
    "cosmic": [r"The evolving flat universe, each point a sphere; $p,q=\arctan[H_0(\eta\mp r)]$. "
               r"The past null boundary has divergent transported tidal curvature [yoshida_quintin2018]."],
    "conformal": [r"The evolving flat universe in conformal time, with $\eta$ extending over the real line. "
                  r"The past null boundary has divergent transported tidal curvature [yoshida_quintin2018]."],
    "closed_de_sitter": [r"The exact closed de Sitter universe, each point a sphere; "
                        r"$T=\arctan[\sinh(H_0\tau)]$, $X=\chi$."],
    "flat_de_sitter": [r"The expanding flat de Sitter patch, with $H_0\eta=-e^{-H_0\tau}$. "
                      r"The past null horizon admits a regular continuation."],
    "open_de_sitter": [r"The expanding open de Sitter patch, with $H_0\eta=\ln\tanh(H_0\tau/2)$. "
                      r"The past light cone admits a regular continuation."],
}


def conformal_diagram(ck, src):
    """Compactify the evolving branch using both sourced asymptotic conformal limits.

    Its past a~exp(tau) and future a~tau^(2/3) imply eta in the whole real line.
    The past edge is explicitly distinguished from complete null infinity.
    For M/H_0=1/2, q<2 and h'/a^2 diverges there; scalar invariants stay finite.
    The exact de Sitter charts use their untruncated conformal-time ranges.
    """
    import conformal as cd
    import slices
    views = []
    pi = np.pi
    for system in CHART_IDS:
        evolving = system in ("cosmic", "conformal")
        radial = r"\chi" if system in ("closed_de_sitter", "open_de_sitter") else "r"
        first = r"\eta" if system == "conformal" else r"\tau"
        params = {"H_0": 1, **({"M": MASS,"k":0} if evolving else {})}
        functions = {"a": "sta_c(eta)" if system == "conformal" else "sta_a(tau)"} if evolving else {}
        plane = cd.Plane(src,ID,system,(first,radial),cd.EQUATOR,params,functions=functions)

        def mapping(t,r,system=system):
            e = conformal_time(system,t)
            return ((e-r)/2,(e+r)/2) if system == "closed_de_sitter" else (np.arctan(e-r),np.arctan(e+r))
        times = ck.uniform(-2,1) if system == "conformal" else ck.uniform(0.3,4) if system == "open_de_sitter" else ck.uniform(-3,3)
        ck.chart(ID+" "+system,plane,mapping,times,ck.uniform(0.02,2.8),lambda t,r:(1,0))
        if system == "closed_de_sitter":
            v = cd.View(system,"Closed de Sitter",[-0.3,pi+0.3,-pi/2-0.3,pi/2+0.3],system)
            square = [[0,-pi/2],[pi,-pi/2],[pi,pi/2],[0,pi/2]]
            v.fill("region",square);v.fill("cover",square)
            v.line("scri",[[[0,-pi/2],[pi,-pi/2]],[[0,pi/2],[pi,pi/2]]])
            v.line("centre",[[[0,-pi/2],[0,pi/2]],[[pi,-pi/2],[pi,pi/2]]])
            for t in (-2,-1,0,1,2):
                v.curve("t",*mapping(np.full(101,t),np.linspace(0,pi,101)))
            v.legend("scri",r"spacelike infinity $\mathscr I^\pm$")
            v.legend("centre","regular poles")
            m = slices.moments(ID,"closed")[0]
            v.slice(m,[mapping(np.zeros(101),np.linspace(0,pi,101))])
        else:
            full = evolving
            bottom,top = -pi, pi if full else 0
            v = cd.View(system,"Evolving flat branch" if full else "Expanding de Sitter",
                        [-0.35,pi+0.35,bottom-0.3,top+0.3],system)
            shape = [[0,bottom],[pi,0],[0,top]]
            v.fill("region",shape);v.fill("cover",shape)
            v.line("centre",[[[0,bottom],[0,top]]])
            v.line("chartedge" if full else "horizon",[[[0,bottom],[pi,0]]])
            v.line("scri",[[[pi,0],[0,top]]])
            v.legend("centre","regular centre")
            v.legend("scri","future null infinity" if full else "future spacelike infinity")
            v.legend("chartedge" if full else "horizon","past null tidal boundary" if full else "regular past horizon")
            evalues = (-4,-2,0,1,4) if full else (-4,-2,-1,-0.3)
            r = cd.spread(0,np.inf,400,10)
            for e in evalues:
                v.curve("t",np.arctan(e-r),np.arctan(e+r))
            if full:
                m = slices.moments(ID,"flat")[0]
                r = np.linspace(0,3,101)
                v.slice(m,[(np.arctan(-r),np.arctan(r))],label=r"$\tau=\eta=0$")
        v.legend("cover","chart domain")
        v.legend("t",r"conformal time constant")
        v.set(settings=r"$H_0=1$"+(r", $M/H_0=1/2$, $k=0$." if evolving else "."))
        if evolving:
            v.set(input=NUMERIC_INPUT)
        views.append(v)
    ck.limit("Starobinsky q<2: the nonzero past mode has divergent h'/a^2",[int(EXPONENT<2)],[1],0.1)
    return views


EMBEDDING_CAPTIONS = {
    "flat": [r"The flat equatorial slice ($\tau=0$, $\theta=\pi/2$), a Euclidean plane with $a=1$."],
    "closed": [r"The closed de Sitter equator ($\tau=0$, $\theta=\pi/2$), a sphere of radius $H_0^{-1}$."],
}


def embedding_diagram(ck, src):
    """One moment of the flat numerical family and the whole exact closed de Sitter equator.

    The flat profile is z=0 and rho=ar.
    The closed profile is rho=sin(chi)/H_0, z=-cos(chi)/H_0 at tau=0.
    Both are read from the published induced metric, through the existing surface integrator.
    """
    import embedding as em
    views = []
    for system,vid in (("cosmic","flat"),("closed_de_sitter","closed")):
        flat = vid == "flat"
        radial = "r" if flat else r"\chi"
        params = {"H_0":1, **({"M":MASS,"k":0} if flat else {})}
        sl = em.Slice(src,ID,system,radial,r"\phi",{"tau":0,"theta":"pi/2"},params,
                      functions={"a":"1"} if flat else None)
        top = 3.0 if flat else np.pi
        size = 6.0 if flat else 2.0
        marks = [(x,"r",None) for x in ((1,2,3) if flat else (np.pi/4,np.pi/2,3*np.pi/4))]
        piece = em.Piece("equator","sheet",sl,0,top,0 if flat else -1,1,
                         (("axis","centre"),("edge","plane continues") if flat else ("axis","antipode")),marks,size)
        surface = em.Surface([piece],label=r"$H_0\tau=0$",time=0.0)
        ck.isometry(ID+" "+vid,piece)
        if flat:
            ck.plane(ID+" flat equator",sl,np.linspace(0.01,3,50))
        else:
            ck.form(ID+" closed sphere",piece,lambda x:-np.cos(x),size)
            ck.radius(ID+" closed circumference",piece,lambda x:np.sin(x),size)
        fig = em.figure_of([surface],{"sheet":"cover"},size,meridians=12)
        fig.legend("fill","cover","equatorial slice")
        fig.legend("line","r",r"comoving circles")
        fig.legend("line","meridian",r"$\phi$ constant")
        views.append(em.view(vid,"Flat slice" if flat else "Closed de Sitter",r"$H_0^{-1}$",
                             [surface],fig.done(),system=system,settings=r"$H_0=1$, $\tau=0$."))
    return views


def curvature_diagram(ck, src):
    """Evaluate the two stored FLRW scalars on the same full-equation trajectory.

    The display includes the de Sitter departure and the first scalaron oscillations.
    No curvature tensor is recomputed, and the axis is ordinary geometric cosmic time.
    """
    import curvature as cv
    params = {"k":0,"H_0":1,"M":MASS}
    functions = {"a":"sta_a(tau)"}
    x_axis = cv.axis(r"$H_0\tau$",-6,80)
    panels,curves = [],[]
    for field,kind,power,lo,hi in (("kretschmann","K",4,0,25),("ricci_scalar","R",2,-1,9)):
        scale = {"kind":"linear","from":lo,"to":hi}
        panels.append({"id":kind,"name":f"${kind}/H_0^{power}$","aspect":2,"scale":scale,
                       "ticks":cv.value_ticks(scale)})
        stored = cv.Published(src,ID,"cosmic",field,["tau"],list(params),functions=functions)
        curves += cv.line_curves(ck,ID+" "+kind,kind,kind,stored,r"\tau",x_axis,scale,cv.PLOT_W/2,params)
    cv.shape_checks(ck,ID,[2,2])
    return [{"id":"evolution","label":"Departure and oscillations","kind":"line","x":x_axis,
             "panels":panels,"curves":curves,"marks":[],
             "legend":[["line","K","$K$"],["line","R","$R$"]],
             "caption":[r"The Kretschmann scalar $K$ and the Ricci scalar $R$ of the evolving flat universe. "
                        r"The Ricci scalar oscillates through zero in the scalaron phase."],
             "ricci_sentences":[[0,1]],"per_scalar":cv.scalar_clips("Linear scale.","Linear scale."),
             "settings":r"$H_0=1$, $M/H_0=1/2$, $k=0$.","input":NUMERIC_INPUT,"clip":"Linear scales."}]
