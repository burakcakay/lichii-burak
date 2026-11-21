"""Deterministic smoke checks; NOT a machine-checked proof of T16/T17."""
import math
import random

def dot(a,b): return sum(t*s for t,s in zip(a,b))
def norm(a): return math.sqrt(dot(a,a))
def vec_add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(s,a):return tuple(s*x for x in a)
def g(a):return scale(norm(a)**-3,a)

r0,R0,b,V=0.30,8.0,0.41,2.0
x=(r0,0.0);u=(0.0,1.5);X=(R0,0.0);W=(3.0,0.0)
d=R0-b
e0=0.5*dot(u,u)-2/r0
ell0=x[0]*u[1]-x[1]*u[0]
s0=dot(X,W)/R0
margins=[b-r0,d,-2/b-e0-2*math.sqrt(b)/(V*d*d),
         abs(ell0)-b*b/(V*d*d),s0*s0-V*V-6/d]
assert all(m>0 for m in margins),margins
lower=(abs(ell0)-b*b/(V*d*d))**2/4
assert 0<lower<r0<b
rng=random.Random(20261121)
for _ in range(500):
    R=rng.uniform(1.5,15.0); bcase=rng.uniform(0.1,min(1.0,R*.7))
    angle=rng.uniform(0,2*math.pi)
    XX=(R*math.cos(angle),R*math.sin(angle))
    rx=rng.uniform(0.0001,bcase)
    phi=rng.uniform(0,2*math.pi)
    xx=(rx*math.cos(phi),rx*math.sin(phi))
    P=g(vec_add(XX,scale(.5,xx)))
    S=g(vec_add(XX,scale(-.5,xx)))
    tide=norm(vec_add(S,scale(-1,P)))
    outer=norm(scale(-1.5,vec_add(P,S)))
    assert tide<=2*rx/(R-bcase)**3*(1+1e-12)
    assert outer<=3/(R-bcase)**2*(1+1e-12)
lam=math.sqrt(6)*(1+math.sqrt(13))/4
ratios=(2*math.sqrt(6)/lam,(math.sqrt(6)/2)/lam,-math.sqrt(6)/lam)
assert 1<ratios[0]<2 and 0<ratios[1]<1 and ratios[2]<0
print("PASS T16/T17 smoke checks; strict margins", [round(v,6) for v in margins],
      "T17 ratios",[round(v,6) for v in ratios])
