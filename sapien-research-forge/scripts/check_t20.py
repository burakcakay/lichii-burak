"""T20 exact local flux tests and T17b/T18a counterexamples.

These tests check finite algebraic identities and numerical accuracy, not
global nonlinear three-body statistics or originality.
"""
import cmath
import json
import math
from pathlib import Path
import random

root=Path(__file__).resolve().parents[1]
packet=json.loads((root/"registry/T17_T20_P20_P27.json").read_text(encoding="utf8"))
assert len(packet["tasks"])==8
assert len(set(x["id"] for x in packet["tasks"]))==8
assert {x["id"] for x in packet["tasks"]}=={f"P{i}" for i in range(20,28)}
assert {x["id"] for x in packet["theorems"]}=={"T20.A","T20.B","T20.C"}

# Check Hermitian rank-one gauge gluing and exact mirror behavior.
rng=random.Random(20261122)
for _ in range(300):
    z=[complex(rng.uniform(-3,3),rng.uniform(-3,3)) for _ in range(4)]
    H=[[z[i].conjugate()*z[j] for j in range(4)] for i in range(4)]
    k=rng.randrange(4);ell=(k+rng.randrange(1,4))%4
    zk=[H[k][j]/math.sqrt(H[k][k].real) for j in range(4)]
    zl=[H[ell][j]/math.sqrt(H[ell][ell].real) for j in range(4)]
    phase=H[k][ell].conjugate()/abs(H[k][ell])
    assert max(abs(phase*zk[j]-zl[j]) for j in range(4)) < 1e-12
    # Flip spatial y signs => z -> conjugate(z), Gram invariant.
    mirror=[t.conjugate() for t in z]
    H_m=[[mirror[i].conjugate()*mirror[j] for j in range(4)] for i in range(4)]
    assert max(abs(H_m[i][j]-H[i][j].conjugate()) for i in range(4) for j in range(4))<1e-12
z=[complex(1,0),complex(0,1)]
h12=z[0].conjugate()*z[1]
zk=[z[0],z[1]]
zl=[-1j,1]
assert all(abs(a-b)<1e-13 for a,b in zip([t*h12.conjugate()/abs(h12) for t in zk],zl))
assert max(abs(a-b) for a,b in zip([t*h12/abs(h12) for t in zk],zl))>1

# Verify exact circular exit flux Jacobian by high-order finite differences.
la,lb,delta,T,psi=3.8,2.1,0.09,2.5,1.2
def point(t,ang):
    return (delta*math.exp(-la*t)*math.cos(ang),
            delta*math.exp(-lb*t)*math.sin(ang))
eps=1e-5
xt=(point(T+eps,psi)[0]-point(T-eps,psi)[0])/(2*eps)
yt=(point(T+eps,psi)[1]-point(T-eps,psi)[1])/(2*eps)
xp=(point(T,psi+eps)[0]-point(T,psi-eps)[0])/(2*eps)
yp=(point(T,psi+eps)[1]-point(T,psi-eps)[1])/(2*eps)
jac=abs(xt*yp-xp*yt)
a=la*math.cos(psi)**2+lb*math.sin(psi)**2
expected=delta**2*math.exp(-(la+lb)*T)*a
assert abs(jac/expected-1)<2e-8,(jac,expected)

# Analytical integral of a(psi): pi*(lambda_plus+lambda_minus).
num=20000
integral=2*math.pi/num*sum(
    la*math.cos(2*math.pi*i/num)**2+
    lb*math.sin(2*math.pi*i/num)**2
    for i in range(num))
assert abs(integral-math.pi*(la+lb))<1e-11

# A smooth positive but anisotropic initial density gives nonuniform finite-R angle law.
R=.1
a_par=.8
p0=(1+2*a_par*R/3)/(2*math.pi)
p_pi=(1-2*a_par*R/3)/(2*math.pi)
assert p0>p_pi>0 and abs((p0+p_pi)/2-1/(2*math.pi))<1e-15

def gamma_and_face_weights(m1,m2,m3):
    sigma=(m1*m2+m1*m3+m2*m3)/(m1+m2+m3)**2
    k=math.sqrt(1-3*sigma)
    Lp,Lm=1+math.sqrt(13+12*k),1+math.sqrt(13-12*k)
    return (Lp+Lm)/4,Lp/(Lp+Lm),Lm/(Lp+Lm)
gamma,wplus,wminus=gamma_and_face_weights(1.,1.,10.)
assert abs(gamma-2.1726039399558577)<1e-12
assert abs(wplus-0.6547921200882851)<1e-12
assert abs(wplus+wminus-1)<1e-15
gamma_eq,wp_eq,wm_eq=gamma_and_face_weights(1.,1.,1.)
assert abs(gamma_eq-(1+math.sqrt(13))/2)<1e-14
assert abs(wp_eq-.5)<1e-14
print("PASS T20 local flux, 300 Hermitian gauge cocycles, finite-R density counterexample")
print("PASS equal-mass gamma",gamma_eq,"unequal gamma",gamma,"exit faces",wplus,wminus)
