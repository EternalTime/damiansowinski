"""Accreting dust shell in an Einstein-de Sitter exterior, G=c=1.

Sakai-Haines (2-8)-(2-10), (2-15), (2-18), (2-20), (2-24)-(2-28).
The angular junction fixes density; the temporal and normal junctions fix motion.
The Newtonian late-time power law is not used as an evolution equation.
"""
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp

T0, R0, V0, END = 1.0, 0.15, 0.0, 4.0


def scale(t):
    return np.asarray(t) ** (2 / 3)


def quantities(t, radius, speed):
    hubble = 2 / (3 * t)
    rho = 1 / (6 * np.pi * t * t)
    gamma = 1 / np.sqrt(1 - speed * speed)
    x = hubble * radius
    clock = np.sqrt(1 + 2 * x * speed + x * x)
    sigma = hubble * hubble * radius / (4 * np.pi * gamma * (clock + 1 + x * speed))
    return hubble, rho, gamma, clock, sigma


def rhs(t, y):
    radius, speed = y[:2]
    hubble, rho, gamma, clock, sigma = quantities(t, radius, speed)
    acceleration = (-gamma * speed * hubble + 2 * np.pi * sigma
                    - rho * gamma * gamma * speed * speed / sigma) / gamma ** 3
    return np.array([hubble * radius + speed, acceleration, clock])


@lru_cache(maxsize=1)
def trajectory():
    solution = solve_ivp(rhs, (T0, END), [R0, V0, 0.0], method="DOP853",
                         rtol=2e-12, atol=2e-14, dense_output=True, max_step=0.025)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def state(t):
    if np.any(np.asarray(t) < T0) or np.any(np.asarray(t) > END):
        raise ValueError("shell trajectory is defined only for 1 <= t <= 4")
    return trajectory().sol(t)



@lru_cache(maxsize=1)
def inner_inverse():
    from scipy.interpolate import PchipInterpolator
    times = np.linspace(T0, END, 8193)
    return PchipInterpolator(state(times)[2], times, extrapolate=False)


def drawn_radius(t):
    times = np.asarray(t)
    valid = (times >= T0) & (times <= END)
    clipped = np.clip(times, T0, END)
    values = state(clipped.ravel())[0].reshape(times.shape)
    return np.where(valid, values, np.nan)


def drawn_inner_radius(T):
    times = inner_inverse()(T)
    return drawn_radius(times)


def declared_functions():
    import sympy as sp
    class sv_radius(sp.Function):
        nargs = 1
        is_real = True
        _imp_ = staticmethod(drawn_radius)
    class sv_inner_radius(sp.Function):
        nargs = 1
        is_real = True
        _imp_ = staticmethod(drawn_inner_radius)
    class sv_comoving_radius(sp.Function):
        nargs = 1
        is_real = True
        _imp_ = staticmethod(lambda t: drawn_radius(t)/scale(t))
    return {"sv_radius": sv_radius, "sv_inner_radius": sv_inner_radius, "sv_comoving_radius": sv_comoving_radius}


DRAWING_INPUT = (r"$t_0=1$, $k=0$, $a=(t/t_0)^{2/3}$; at $t=t_0$, $R=0.15t_0$, $v=0$, and $T=0$.")
DRAWING_FUNCTIONS = {"a": "t**Rational(2,3)", "R": "sv_radius(t)"}


def diagram_specs(Diagram, EQUATOR, SIDEWAYS):
    return [
        Diagram("shell_void", "interior", "through", "Through the cavity", ("T", "r"),
                (0, .95, 0, float(state(END)[2])), "$x/t_0$", "$T/t_0$", {}, EQUATOR,
                mirror=True, tau="T", families=SIDEWAYS, cones=(3,9), areal=True,
                functions={"R": "sv_inner_radius(T)"}, where="sv_inner_radius(T)-r", input=DRAWING_INPUT,
                curves=(("shell", "sv_inner_radius(T)", "the shell"),)),
        Diagram("shell_void", "exterior", "radial", "Outside the shell", ("t", "r"),
                (0, 1.9, T0, END), "$r/t_0$", "$t/t_0$", {"k": 0}, EQUATOR,
                tau="t", cones=(6,9), areal=True, functions=DRAWING_FUNCTIONS,
                where="r-sv_radius(t)/t**Rational(2,3)", input=DRAWING_INPUT,
                curves=(("shell", "sv_radius(t)/t**Rational(2,3)", "the shell"),)),
        Diagram("shell_void", "exterior_proper", "radial", "Comoving radial coordinate", ("t", r"\chi"),
                (0, 1.9, T0, END), r"$\chi$", "$t/t_0$", {"k":0}, EQUATOR,
                tau="t", cones=(6,9), areal=True,
                functions={"a":"t**Rational(2,3)", "f":"chi", "chi_s":"sv_comoving_radius(t)"},
                where="chi-sv_comoving_radius(t)",
                input=DRAWING_INPUT.replace("$a=(t/t_0)^{2/3}$", "$a=t_0(t/t_0)^{2/3}$"),
                curves=(("shell","sv_comoving_radius(t)","the shell"),)),
    ]


def diagram_captions():
    return {
        ("shell_void", "interior", "through"): [r"The radial plane through the empty cavity ($x=\pm r$, $\theta=\pi/2$, $\phi=0$), bounded by the moving shell."],
        ("shell_void", "exterior_proper", "radial"): [r"The radial plane ($\theta=\pi/2$, $\phi=0$) in the flat dust exterior, with horizontal distance $t_0\chi$."],
        ("shell_void", "exterior", "radial"): [r"The radial plane of the dust exterior ($\theta=\pi/2$, $\phi=0$), with the shell sweeping outward in comoving radius."],
    }


def null_forms():
    """Dimensionless cavity null coordinates in units of its initial radius.

    The cavity diameter stays below one t0 in the interval drawn; using R0 as
    the null-coordinate unit makes the family separation test measure its size.
    """
    return {
        ("shell_void", "interior", "through"): (lambda T,r:(T+r)/R0, lambda T,r:(T-r)/R0, None),
        ("shell_void", "exterior_proper", "radial"): (lambda t,r:3*np.cbrt(t)+r, lambda t,r:3*np.cbrt(t)-r, None),
        ("shell_void", "exterior", "radial"): (lambda t,r:3*np.cbrt(t)+r, lambda t,r:3*np.cbrt(t)-r, None),
    }


def embedding_views(ck, src, api):
    """Piecewise constant T and t at one shell event, intrinsically Euclidean on each side.

    The spacelike slice has a corner in its time normal at the shell; its intrinsic
    radial distance and circumference agree there. It has no intrinsic throat.
    """
    time = 2.0
    R,v,T = state(time)
    a = float(scale(time))
    top, size = .75, 1.5
    flat = api.Slice(src,"shell_void","interior","r",r"\phi",{"T":float(T),**api.EQUATOR},{})
    outer = api.Slice(src,"shell_void","exterior","r",r"\phi",{"t":time,**api.EQUATOR},
                      {"k":0}, {"a":"2**Rational(2,3)"})
    inside = api.Piece("inside","sheet",flat,0.0,float(R),0.0,1,
                       (("axis",r"the centre $r=0$"),("join","the shell")),[(float(R),"surface",None)],size)
    outside = api.Piece("outside","sheet",outer,float(R/a),top/a,0.0,1,
                        (("join","the shell"),("edge","the slice continues outward")),[],size)
    ck.plane("shell void, the cavity slice",flat,np.linspace(.001,R,100))
    ck.plane("shell void, the exterior slice",outer,np.linspace(R/a,top/a,100))
    for piece in (inside,outside):
        ck.isometry("shell void, "+piece.id,piece)
    ck.join("shell void, the two pieces meet",inside,float(R),outside,float(R/a))
    surface = api.Surface([inside,outside],label="$t=2t_0$",time=time)
    fig = api.figure_of([surface],{"sheet":"cover"},size)
    api.ring_label(fig,[0,0,0],float(R),0,r"$R=0.244t_0$",dx=12)
    fig.legend("fill","cover","the spatial slice")
    fig.legend("line","surface","the shell")
    fig.legend("line","meridian",r"$\phi$ constant")
    figure = fig.done()
    x0,x1,y0,y1 = figure["box"]
    half = max(x1-x0,y1-y0)/2
    centre = (y0+y1)/2
    figure["box"] = [x0,x1,centre-half,centre+half]
    return [api.view("slice","The cavity and shell","$t_0$",[surface],figure,
                     settings=r"$t_0=1$; $t=2t_0$ outside, $T=1.00439t_0$ inside.",input=DRAWING_INPUT)]


def embedding_captions():
    return {("shell_void","slice"):[r"The flat equatorial slice ($\theta=\pi/2$), with the shell at $R=0.244t_0$."]}


def flat_slice_marks(api,system):
    out=[]
    for moment in api.moments("shell_void"):
        lo,hi=moment.reach("exterior" if system=="exterior_proper" else system,"r")
        t=float(state(2.0)[2]) if system=="interior" else 2.0
        out.append(api.Mark(moment,[[(t,lo),(t,hi)]]))
    return out


def curvature_views(ck,src,api):
    """Stored smooth exterior K only, with a blank cavity and the shell marked separately.

    No product of surface delta distributions is formed. The cosmological scale
    factor specializes the already stored scalar; curvature is not recomputed.
    """
    x,y=api.axis("$r/t_0$",0,1.9),api.axis("$t/t_0$",T0,END)
    K=api.nr_published(src,"shell_void","exterior","kretschmann",["r","t"],["k"],functions=DRAWING_FUNCTIONS)
    nx,ny=95,150
    xs,ys,grid=api.cells(K,x,y,nx,ny,[0])
    X,Y=np.meshgrid(xs,ys)
    grid=np.broadcast_to(grid,X.shape)
    allowed=X>drawn_radius(Y)/scale(Y)
    rows=[[api.rounded(value) if ok else None for value,ok in zip(row,mask)] for row,mask in zip(grid,allowed)]
    times=np.linspace(T0,END,201)
    world=[[float(R/scale(t)),float(t)] for t,R in zip(times,drawn_radius(times))]
    return [{"id":"bulk","label":"Smooth exterior curvature","system":"exterior","kind":"map","x":x,"y":y,
             "colour":{"name":"$Kt_0^4$","scale":{"kind":"log","from":.01,"to":3.0},
                       "ticks":[{"at":.01,"label":"$0.01$"},{"at":.1,"label":"$0.1$"},{"at":1.,"label":"$1$"}]},
             "map":{"metric":"shell_void","system":"exterior","scalar":"kretschmann","x":"r","y":"t",
                    "params":{"k":0},"functions":dict(DRAWING_FUNCTIONS),
                    "domain":"r>sv_radius(t)/t**Rational(2,3)","nx":nx,"ny":ny,"cells":rows},
             "marks":[{"class":"front","points":world}],"legend":[["line","front","the shell"]],
             "caption":[r"The smooth exterior Kretschmann scalar $K$, with the flat cavity blank and the distributional shell marked by its trajectory."],
             "settings":DRAWING_INPUT,"clip":"Logarithmic scale."}]
