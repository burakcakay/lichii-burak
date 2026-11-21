"""T18: finite conservative X10 box-containment smoke and random sampling."""
import math
import random

def vlen(v):return math.hypot(*v)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def margin(z,b,V):
    x,X,u,W=z
    r=vlen(x);R=vlen(X);d=R-b
    e=dot(u,u)/2-2/r
    ell=abs(x[0]*u[1]-x[1]*u[0])
    s=dot(X,W)/R
    return [b-r,d,-2/b-e-2*math.sqrt(b)/(V*d*d),
            ell-b*b/(V*d*d),s*s-V*V-6/d,s]
def box_margins(x,X,u,W,eps,b,V):
    delta=math.sqrt(2)*eps
    rb=vlen(x);RB=vlen(X);ub=vlen(u);Wb=vlen(W)
    rplus=rb+delta;Rminus=RB-delta;Rplus=RB+delta
    dminus=Rminus-b
    eplus=.5*(ub+delta)**2-2/rplus
    Lminus=abs(x[0]*u[1]-x[1]*u[0])-delta*(rb+ub)-delta**2
    Nminus=dot(X,W)-delta*(RB+Wb)-delta**2
    sminus=Nminus/Rplus
    return [rb-delta,b-rplus,dminus,Nminus,-2/b-eplus-2*math.sqrt(b)/(V*dminus*dminus),
            Lminus-b*b/(V*dminus*dminus),sminus*sminus-V*V-6/dminus]
x=(.3,0.);X=(8.,0.);u=(0.,1.5);W=(3.,0.)
b,V,eps=.41,2.,.02
M=box_margins(x,X,u,W,eps,b,V)
assert min(M)>0.035
rng=random.Random(20261121)
for _ in range(3000):
    point=tuple(tuple(val+rng.uniform(-eps,eps) for val in block) for block in (x,X,u,W))
    assert min(margin(point,b,V))>0
print("PASS T18 conservative full 8D cube certificate; min margin",min(M))
