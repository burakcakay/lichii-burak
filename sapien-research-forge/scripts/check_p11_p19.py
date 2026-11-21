"""Check P11-P19 research queue invariants; not a machine proof of conjectures."""
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
P = json.loads((ROOT/"registry/G03_T16_P11_P19.json").read_text(encoding="utf-8"))
assert P["artifact_id"] == "THREEBODY_G03_T16_BRIDGE_P11_P19_V1"
nodes = {item["id"]: item for item in P["problems"]}
assert set(nodes) == {f"P{i}" for i in range(11,20)}
assert len(P["critical_path"]) == len(set(P["critical_path"]))
for n in nodes.values():
    assert n["unknowns"] and n["equations"] and n["deliverable"]
    assert n["accept"] and n["falsify"]
    assert set(n["depends_on"]) <= set(nodes)

visiting=set();visited=set()
def visit(node):
    if node in visiting: raise ValueError("cyclic dependency: " + node)
    if node in visited: return
    visiting.add(node)
    for dep in nodes[node]["depends_on"]:visit(dep)
    visiting.remove(node);visited.add(node)
for node in nodes:visit(node)

s = (1064 - 32*math.sqrt(129))/5625
k = (48 + math.sqrt(129))/75
a, b = math.sqrt(13+12*k),math.sqrt(13-12*k)
assert abs(k*k -(1-3*s)) < 1e-14
assert abs(a-2*b-1) < 1e-13
lambdas=np.array([1+a,1+b,1-a,1-b])
assert abs(lambdas[1]+lambdas[2]-lambdas[3]) < 1e-13

x = np.array([0.3,0.0]);u=np.array([0.0,1.5])
X = np.array([8.0,0.0]);W=np.array([3.0,-27/640])
total_L = .5*np.linalg.det(np.stack([x,u])) + (2/3)*np.linalg.det(np.stack([X,W]))
assert abs(total_L)<1e-14
v=np.array([x/math.sqrt(2),math.sqrt(2/3)*X,u/math.sqrt(2),math.sqrt(2/3)*W])
J2=np.array([[0.,1.],[-1.,0.]])
Gamma=v@v.T; Omega=v@J2@v.T
v_mirror=v@np.diag([1.,-1.])
assert np.max(abs(Gamma-v_mirror@v_mirror.T))<1e-14
assert np.max(abs(Omega+v_mirror@J2@v_mirror.T))<1e-14
assert abs(Omega[0,2]+Omega[1,3]) < 1e-14
bb=.41;V=2.;d=np.linalg.norm(X)-bb
e0=.5*(u@u)-2/np.linalg.norm(x)
ell0=np.linalg.det(np.stack([x,u]))
ss=X@W/np.linalg.norm(X)
margins=[bb-np.linalg.norm(x),d,-2/bb-e0-2*math.sqrt(bb)/(V*d*d),
         abs(ell0)-bb*bb/(V*d*d),ss*ss-V*V-6/d]
assert min(margins)>0.1
assert ell0-bb*bb/(V*d*d)>0

def scaled(z):
    x,X,u,W=z[:2],z[2:4],z[4:6],z[6:8]
    nu=.5*x@u+(2/3)*X@W
    pp=X+.5*x;ss=X-.5*x
    P=pp/np.linalg.norm(pp)**3
    S=ss/np.linalg.norm(ss)**3
    return np.r_[u-nu*x,W-nu*X,
                 -2*x/np.linalg.norm(x)**3+S-P+.5*nu*u,
                 -1.5*(P+S)+.5*nu*W]

nuE=math.sqrt(5*math.sqrt(2))
A=math.sqrt(237/5)
B=math.sqrt(51/5)
zE=np.r_[math.sqrt(2),0.,0.,0.,nuE*math.sqrt(2),0.,0.,0.]
assert np.max(abs(scaled(zE)))<1e-13
h=5e-7
JJ=np.column_stack([(scaled(zE+h*np.eye(8)[i])-scaled(zE-h*np.eye(8)[i]))/(2*h)
                    for i in range(8)])
obs=np.linalg.eigvals(JJ)
expected=np.array([0,nuE,-nuE/2,-2*nuE,
                   nuE*(A-1)/4,-nuE*(A+1)/4,
                   nuE*(-1+1j*B)/4,nuE*(-1-1j*B)/4],dtype=complex)
row,col=linear_sum_assignment(abs(obs[:,None]-expected[None,:]))
err=max(abs(obs[row]-expected[col]))
assert err<2e-7,err
print("PASS P11-P19 schema/DAG, sigma* exact resonance, zero-L mirror escape and Euler eigen test")
print("Euler nuE",nuE,"spectrum residual",err,"X10 min margin",min(margins))
