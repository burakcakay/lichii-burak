#!/usr/bin/env python3
"""Uc cisim: Tur 3 bagimsiz denetim.

S1: E>0 iken global D=0 kesitine karsi acik analitik ornek.
S2: Kütleli Wronskiyen sabitlerinin sayisal denetimi.
S3: Gercekten 13 skalerli polinomsal, rotasyondan arindirilmis dinamik.
S4: Uclu carpisma yakininda olcek uslerinin analitik/sayisal karsilastirmasi.

Gereksinimler: numpy, scipy. Calistir: python uc_cisim_tur3_dogrulama.py
S4'u atlamak icin: python uc_cisim_tur3_dogrulama.py --no-s4
Sayisal sonuclar ispat degildir; S4'te hata kontrollu keyfi duyarlik yoktur.
"""
import argparse
import math
import numpy as np
from scipy.integrate import solve_ivp

PAIRS=((0,1),(0,2),(1,2))
RTOL=1e-12
ATOL=1e-12

def newton_acc(r,m,G=1.):
    a=np.zeros_like(r)
    for i,j in PAIRS:
        d=r[j]-r[i]
        f=G*d/np.linalg.norm(d)**3
        a[i]+=m[j]*f
        a[j]-=m[i]*f
    return a

def newton(t,z,m,G=1.):
    r=z[:9].reshape(3,3);v=z[9:].reshape(3,3)
    return np.r_[v.ravel(),newton_acc(r,m,G).ravel()]

def xrel(r):
    return np.stack([r[2]-r[1],r[0]-r[2],r[1]-r[0]],axis=1)

def pack13(r,v):
    B=xrel(r)[:,:2]; W=xrel(v)[:,:2]
    Q=B.T@B;D=B.T@W;P=W.T@W
    u=1/np.linalg.norm(xrel(r),axis=0)
    return np.r_[Q[0,0],Q[0,1],Q[1,1],D.ravel(),P[0,0],P[0,1],P[1,1],u]

def unpack13(y):
    Q=np.array([[y[0],y[1]],[y[1],y[2]]]);D=y[3:7].reshape(2,2)
    P=np.array([[y[7],y[8]],[y[8],y[9]]]);u=y[10:13]
    return Q,D,P,u

def reduction13(t,y,m,G=1.):
    Q,D,P,u=unpack13(y)
    m1,m2,m3=m
    w=G*u**3
    A=np.array([
        [-(m2+m3)*w[0]-m1*w[2], m2*(w[0]-w[2])],
        [m1*(w[1]-w[2]), -(m1+m3)*w[1]-m2*w[2]]
    ])
    Qd=D+D.T
    Dd=P+Q@A
    Pd=A.T@D+D.T@A
    ud=-u**3*np.array([D[0,0],D[1,1],D.sum()])
    return np.r_[Qd[0,0],Qd[0,1],Qd[1,1],Dd.ravel(),
                 Pd[0,0],Pd[0,1],Pd[1,1],ud]

def invariant13(y,m,G=1.):
    Q,D,P,u=unpack13(y)
    mu=np.array([m[1]*m[2],m[0]*m[2],m[0]*m[1]])
    M=sum(m)
    H=np.array([[mu[0]+mu[2],mu[2]],[mu[2],mu[1]+mu[2]]])
    energy=np.sum(H*P)/(2*M)-G*np.dot(mu,u)
    L2=sum(H[a,b]*H[c,d]*(Q[a,c]*P[b,d]-D[a,d]*D[c,b])
           for a in range(2) for b in range(2)
           for c in range(2) for d in range(2))/M**2
    return energy,L2

def test_s1():
    # Side=1; mass=1; homotetik hareket r_i=a(t)r_i0,
    # a''=-3/a^2. a(0)=1, a'(0)=3.
    # E=(1/2)a'^2 - 3/a = 3/2; D=a*a'>0 until past triple collision.
    E=.5*3**2-3
    print('S1: global section counterexample: E=',E,
          'D(0)=',3., 'regular backward D=0 intersection = NONE')

def test_s2():
    m=np.array([1.,2.5,.7]);rng=np.random.default_rng(3)
    r0=rng.normal(size=(3,3));v0=.35*rng.normal(size=(3,3))
    q0=np.random.default_rng(11).normal(size=3)
    qd0=np.random.default_rng(12).normal(size=3)
    def coupled(t,y):
        r=y[:9].reshape(3,3);v=y[9:18].reshape(3,3)
        q=y[18:21];qd=y[21:24];qdd=np.zeros(3)
        for i,j in PAIRS:
            s=(q[j]-q[i])/np.linalg.norm(r[j]-r[i])**3
            qdd[i]+=m[j]*s;qdd[j]-=m[i]*s
        return np.r_[v.ravel(),newton_acc(r,m).ravel(),qd,qdd]
    sol=solve_ivp(coupled,[0,4],np.r_[r0.ravel(),v0.ravel(),q0,qd0],
                  method='DOP853',rtol=RTOL,atol=ATOL,dense_output=True)
    assert sol.success
    def wronsk(t):
        y=sol.sol(t);r=y[:9].reshape(3,3);v=y[9:18].reshape(3,3)
        q=y[18:21];qd=y[21:24]
        basis=[(r[:,k],v[:,k]) for k in range(3)]+[(np.ones(3),np.zeros(3)),(t*np.ones(3),np.ones(3))]
        W=lambda a,ad,b,bd:np.sum(m*(a*bd-ad*b))
        return np.r_[[W(q,qd,b,bd) for b,bd in basis],
                     [W(*basis[i],*basis[j]) for i in range(5) for j in range(i+1,5)]]
    start=wronsk(0)
    drift=max(np.max(np.abs(wronsk(t)-start)) for t in np.linspace(0,4,41))
    print('S2: max Wronskian drift =',f'{drift:.3e}')

def test_s3():
    m=np.array([1.,2.5,.7]); rng=np.random.default_rng(3)
    r0=rng.normal(size=(3,3));v0=.35*rng.normal(size=(3,3))
    y0=pack13(r0,v0)
    full=solve_ivp(newton,[0,4],np.r_[r0.ravel(),v0.ravel()],args=(m,),
                   method='DOP853',rtol=RTOL,atol=ATOL)
    red=solve_ivp(reduction13,[0,4],y0,args=(m,),
                  method='DOP853',rtol=RTOL,atol=ATOL)
    assert full.success and red.success
    z=full.y[:,-1]
    target=pack13(z[:9].reshape(3,3),z[9:].reshape(3,3))
    actual=red.y[:,-1]
    q,d,p,u=unpack13(actual)
    check=np.array([u[0]**2*q[0,0]-1,u[1]**2*q[1,1]-1,
                    u[2]**2*(q[0,0]+2*q[0,1]+q[1,1])-1])
    drift=np.array(invariant13(actual,m))-np.array(invariant13(y0,m))
    gram_det=np.linalg.det(np.block([[q,d],[d.T,p]]))
    print('S3: scalar ODE states =',len(y0),
          '; max |13-Newton| =',f'{np.max(np.abs(actual-target)):.3e}')
    print('S3: invariant drift E,L2=',drift,
          '; max norm constraint=',f'{max(np.abs(check)):.3e}',
          '; Gram det=',f'{gram_det:.3e}')

def s4_rhs(t,y):
    r=y[:6].reshape(3,2);v=y[6:].reshape(3,2)
    a=np.zeros((3,2))
    for i,j in PAIRS:
        d=r[j]-r[i]
        f=d/np.linalg.norm(d)**3
        a[i]+=f;a[j]-=f
    return np.r_[v.ravel(),a.ravel()]

def collapse_sample(eps,phi):
    ang=np.array([np.pi/2,7*np.pi/6,11*np.pi/6]);R=1/math.sqrt(3)
    r=R*np.stack([np.cos(ang),np.sin(ang)],axis=1)
    theta=np.pi/2+np.deg2rad(phi)
    r[0]+=eps*np.array([math.cos(theta),math.sin(theta)])
    r-=r.mean(axis=0)
    z0=np.r_[r.ravel(),np.zeros(6)]
    def turn(t,y):return np.dot(y[:6],y[6:])
    turn.terminal=True;turn.direction=1
    kw=dict(method='DOP853',rtol=2e-12,atol=2e-12,max_step=.01)
    sol=solve_ivp(s4_rhs,[0,3],z0,events=turn,**kw)
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('No turnaround')
    z=sol.y_events[0][0]
    rmin=np.linalg.norm(z[:6])
    def stop(t,y):return np.linalg.norm(y[:6])-10*rmin
    stop.terminal=True;stop.direction=1
    end=solve_ivp(s4_rhs,[sol.t_events[0][0],3],z,events=stop,**kw)
    if not end.success or not len(end.t_events[0]):raise RuntimeError('No exit')
    ze=end.y_events[0][0]; rr=ze[:6].reshape(3,2);vv=ze[6:].reshape(3,2)
    dd=np.array([np.linalg.norm(rr[j]-rr[i]) for i,j in PAIRS])
    i,j=PAIRS[np.argmin(dd)];k=3-i-j
    relv=vv[k]-(vv[i]+vv[j])/2
    return rmin,np.linalg.norm(relv),k+1

def test_s4():
    e_r=4/(1+math.sqrt(13))
    e_v=2/(1+math.sqrt(13))
    print('S4: theoretical exponents: radius=',f'{e_r:.9f}',
          ', velocity=',f'{-e_v:.9f}')
    for phi in (40,100):
        rows=[]
        for eps in (1e-2,1e-3,1e-4):
            radius,speed,label=collapse_sample(eps,phi)
            rows.append((radius,speed))
            print(f'S4 phi={phi:3d}, eps={eps:.0e}: rmin={radius:.9e}, '
                  f'v={speed:.9f}, v*sqrt(rmin)={speed*math.sqrt(radius):.9f}, '
                  f'candidate={label}')
        r_exp=math.log10(rows[1][0]/rows[2][0])
        v_exp=math.log10(rows[2][1]/rows[1][1])
        print(f'S4 phi={phi}: last-decade measured (radius, velocity) exponents = '
              f'({r_exp:.6f}, {-v_exp:.6f})')


# 9 skalere ek yerel indirgeme (kollinear-olmayan durumlar):
# K = P - D^T Q^{-1} D = h h^T, K rank <= 1,
# h' = -D^T Q^{-1} h; P bundan geri kurulur.
def pack9(r, v):
    B=xrel(r)[:,:2];W=xrel(v)[:,:2]
    Q=B.T@B;D=B.T@W;P=W.T@W
    K=(P-D.T@np.linalg.solve(Q,D))
    K=(K+K.T)/2
    e,U=np.linalg.eigh(K)
    h=np.sqrt(max(e[-1],0))*U[:,-1]
    return np.r_[Q[0,0],Q[0,1],Q[1,1],D.ravel(),h]

def reduction9(t, y, m, G=1.):
    Q=np.array([[y[0],y[1]],[y[1],y[2]]]);D=y[3:7].reshape(2,2);h=y[7:9]
    if np.linalg.det(Q)<=0: raise ValueError('Collinear/singular shape, switch chart')
    m1,m2,m3=m
    u=np.array([1/np.sqrt(Q[0,0]),1/np.sqrt(Q[1,1]),1/np.sqrt(Q.sum())])
    w=G*u**3
    A=np.array([[ -(m2+m3)*w[0]-m1*w[2], m2*(w[0]-w[2]) ],
                [ m1*(w[1]-w[2]), -(m1+m3)*w[1]-m2*w[2] ]])
    QinvD=np.linalg.solve(Q,D)
    P=D.T@QinvD+np.outer(h,h)
    qdot=D+D.T
    ddot=P+Q@A
    hdot=-QinvD.T@h
    return np.r_[qdot[0,0],qdot[0,1],qdot[1,1],ddot.ravel(),hdot]

def test_s3_nine():
    m=np.array([1.,2.5,.7]);rng=np.random.default_rng(3)
    r0=rng.normal(size=(3,3));v0=.35*rng.normal(size=(3,3))
    first=pack9(r0,v0)
    a=solve_ivp(reduction9,[0,4],first,args=(m,),method='DOP853',rtol=1e-13,atol=1e-13)
    b=solve_ivp(newton,[0,4],np.r_[r0.ravel(),v0.ravel()],args=(m,),method='DOP853',rtol=1e-13,atol=1e-13)
    assert a.success and b.success
    b9=pack9(b.y[:9,-1].reshape(3,3),b.y[9:,-1].reshape(3,3))
    # The square root factor h is unique only up to a sign.
    direct=np.max(np.abs(a.y[:7,-1]-b9[:7]))
    herror=min(np.linalg.norm(a.y[7:,-1]-b9[7:]),np.linalg.norm(a.y[7:,-1]+b9[7:]))
    print('S3-9:', 'Q,D max abs mismatch',direct,'h +/- error',herror,
          'minimum det Q',min(np.linalg.det(np.array([[yy[0],yy[1]],[yy[1],yy[2]]])) for yy in a.y.T))

def reconstruct13free(y9, E, m, Rcm=None, Vcm=None):
    B=E[:,:2];n=E[:,2]
    Q=np.array([[y9[0],y9[1]],[y9[1],y9[2]]])
    D=y9[3:7].reshape(2,2);h=y9[7:9]
    Bd=B@np.linalg.solve(Q,D)+np.outer(n,h)
    x1,x2=B[:,0],B[:,1]
    w1,w2=Bd[:,0],Bd[:,1]
    M=sum(m)
    r2=-((m[0]+m[2])*x1+m[0]*x2)/M
    v2=-((m[0]+m[2])*w1+m[0]*w2)/M
    r=np.stack([r2+x1+x2,r2,r2+x1],axis=0)
    v=np.stack([v2+w1+w2,v2,v2+w1],axis=0)
    if Rcm is not None:r+=Rcm
    if Vcm is not None:v+=Vcm
    return r,v

def test_full_reconstruction():
    m=np.array([1.,2.5,.7]);rng=np.random.default_rng(3)
    r0=rng.normal(size=(3,3));v0=.35*rng.normal(size=(3,3))
    B=xrel(r0)[:,:2];Bd=xrel(v0)[:,:2]
    n=np.cross(B[:,0],B[:,1]);n/=np.linalg.norm(n)
    E0=np.column_stack((B,n));h=Bd.T@n
    y0=pack9(r0,v0);y0[7:9]=h
    R0=np.average(r0,axis=0,weights=m);V0=np.average(v0,axis=0,weights=m)
    def lifted(t,Y):
        y=Y[:9];E=Y[9:].reshape(3,3)
        Q=np.array([[y[0],y[1]],[y[1],y[2]]]);D=y[3:7].reshape(2,2);h=y[7:9]
        qd=np.linalg.solve(Q,D)
        W=np.block([[qd, -np.linalg.solve(Q,h).reshape(2,1)],
                    [h.reshape(1,2),np.zeros((1,1))]])
        return np.r_[reduction9(t,y,m), (E@W).ravel()]
    lifted_sol=solve_ivp(lifted,[0,4],np.r_[y0,E0.ravel()],method='DOP853',rtol=1e-13,atol=1e-13)
    direct=solve_ivp(newton,[0,4],np.r_[r0.ravel(),v0.ravel()],args=(m,),method='DOP853',rtol=1e-13,atol=1e-13)
    Y=lifted_sol.y[:,-1]
    r,v=reconstruct13free(Y[:9],Y[9:].reshape(3,3),m,R0+4*V0,V0)
    expected=direct.y[:,-1]
    print('9D shape + linear frame reconstruction: max |r-rN|',np.max(abs(r.ravel()-expected[:9])),
          'max |v-vN|',np.max(abs(v.ravel()-expected[9:])),
          'frame normal error',abs(np.linalg.norm(Y[9:].reshape(3,3)[:,2])-1))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--no-s4',action='store_true')
    args=parser.parse_args()
    test_s1()
    test_s2()
    test_s3()
    test_s3_nine()
    test_full_reconstruction()
    if not args.no_s4:
        test_s4()