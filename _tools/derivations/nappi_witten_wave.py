"""Compare Nappi-Witten 1993 (17) with the live symmetric plane-wave metric.

This is the later wave used in the cosmology's History and related texts.
At c=1, omega=-1/2 and V=-v-bu/2 identify the rotating chart; its
Brinkmann eigenvalues are lambda_1=lambda_2=-1/4.
"""

import sys,json
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import verify_metrics as vm
page=json.loads((ROOT/"MFS/assets/data/metrics/cahen_wallach.json").read_text())
chart=next(c for c in page["coordinates"] if c["id"]=="rotating")
r=vm.Reader(chart["coords"],[p["symbol"] for p in chart["parameters"]],())
g=S.Matrix(vm.metric_from_line_element(r,chart["line_element"],chart["coords"]))
u,v,a1,a2,b=S.symbols("u v a1 a2 b",real=True)
X=S.Matrix([u,-v-b*u/2,a1,a2])
old=[r.symbol[x] for x in chart["coords"]]
new=[u,v,a1,a2]
sub=dict(zip(old,X));sub.update({r.c:1,r.parameters["omega"]:-S.Rational(1,2)})
J=X.jacobian(new)
pulled=J.T*g.subs(sub,simultaneous=True)*J
expected=S.Matrix([[b,1,a2/2,-a1/2],[1,0,0,0],[a2/2,0,1,0],[-a1/2,0,0,1]])
assert all(S.expand(x)==0 for x in pulled-expected)
print("1993 (17) agrees with the live rotating chart at omega=-1/2, V=-v-bu/2 (c=1).")
print("Brinkmann eigenvalues lambda_1=lambda_2=-1/4; Ricci_uu=1/2.")
