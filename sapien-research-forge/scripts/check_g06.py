"""G06: exact source-density and endpoint toy-model identities; NOT a three-body CAP."""
from math import sqrt,log
from itertools import combinations_with_replacement
import sympy as sp

rates=((1+sqrt(22))/4,0.75)
gamma=sum(rates)
assert abs(gamma-2.1726039399558577)<1e-14
biased_gamma=3*rates[0]+rates[1]
assert abs(biased_gamma-5.017811819867572)<1e-14
assert abs(3*rates[0]/biased_gamma-0.8505324577875872)<1e-14

def tail_integral(r,g,p):
    a=r**(1/p)
    if abs(p*g-1)<1e-12:return r**g*(1+log(1/r)/p)
    return a+r**g*(1-a**(1-p*g))/(1-p*g)
for p,expected in [(0.2,gamma),(1.0,1.0)]:
    a,b=1e-15,1e-18
    measured=log(tail_integral(b,gamma,p)/tail_integral(a,gamma,p))/log(b/a)
    assert abs(measured-expected)<0.02,(p,measured,expected)

# Exhaustive exact 40 quadratic eigenvalue checks, on physical 0<k<1.
cs=[(1,0),(0,1),(-1,0),(0,-1)]
relations=[]
for i in range(4):
    for j,l in combinations_with_replacement(range(4),2):
        P=cs[i][0]-cs[j][0]-cs[l][0]
        Q=cs[i][1]-cs[j][1]-cs[l][1]
        if P==0: continue  # Q*b=1 has no solution when b>1.
        for eps in (-1,1):
            b=(sp.Integer(Q)+eps*abs(P)*sp.sqrt(26*(P*P+Q*Q)-1))/(P*P+Q*Q)
            a=(1-Q*b)/P
            if 1<float(b)<sqrt(13) and sqrt(13)<float(a)<5:
                k=sp.simplify((13-b*b)/12)
                relations.append((i+1,j+1,l+1,k))
kp=(48+sp.sqrt(129))/75
km=(48-sp.sqrt(129))/75
assert len(relations)==4
assert all(k in (kp,km) for _,_,_,k in relations)
assert {(i,j,l) for i,j,l,k in relations if k==kp}=={(1,2,2),(4,2,3)}
assert {(i,j,l) for i,j,l,k in relations if k==km}=={(2,1,4),(3,4,4)}
print("PASS G06 toy exact tails and 40-case exact quadratic resonance enumeration")
