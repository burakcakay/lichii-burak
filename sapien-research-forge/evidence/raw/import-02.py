#!/usr/bin/env python3
"""Uc cisim, tur 4: Gram matrisinde yapi-koruyan donmus-katsayi ustelli integrator.

S'=K(w)^T S + S K(w), K = [[0,A(w)],[I,0]].
Her sabit w adiminda S_new=exp(hK).T @ S_old @ exp(hK):
- det S, rank S ve S>=0 yapisi cebirsel olarak korunur
- uzaysal acisal momentum |L|^2 korunur, cunku A H = H A^T
- enerji ise YALNIZCA yaklasik korunur
Sinif: 13-degiskenli surekli sistemin 10-boyutlu Gram alt yapisini,
mesafe kisitlarini cebirsel olarak geri hesaplayarak entegre eder.
Gereksinimler: numpy, scipy.
"""
import argparse
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp


def relative(r):
    return np.stack((r[2]-r[1], r[0]-r[2]), axis=1)


def Hmat(m):
    m1,m2,m3=m
    return np.array([[m2*(m1+m3),m1*m2],[m1*m2,m1*(m2+m3)]],float)


def Amat(m, w):
    m1,m2,m3=m
    w1,w2,w3=w
    return np.array([[-(m2+m3)*w1-m1*w3, m2*(w1-w3)],
                     [m1*(w2-w3), -(m1+m3)*w2-m2*w3]],float)


def Kmat(m, w):
    K=np.zeros((4,4),float)
    K[:2,2:]=Amat(m,w)
    K[2:,:2]=np.eye(2)
    return K


def getw(S,G=1.0):
    Q=S[:2,:2]
    d2=np.array([Q[0,0],Q[1,1],Q[0,0]+2*Q[0,1]+Q[1,1]])
    if min(d2)<=0:
        raise FloatingPointError('Kare mesafeler pozitif degil')
    return G/d2**1.5


def kinetic_potential(S,m,G=1.0):
    H=Hmat(m);M=m.sum()
    T=np.sum(H*S[2:,2:])/(2*M)
    mu=np.array([m[1]*m[2],m[0]*m[2],m[0]*m[1]])
    Q=S[:2,:2]
    d=np.sqrt(np.array([Q[0,0],Q[1,1],Q[0,0]+2*Q[0,1]+Q[1,1]]))
    U=G*np.dot(mu,1/d)
    return T-U


def L2_from_gram(S,m):
    H=Hmat(m);M=m.sum()
    Q=S[:2,:2];D=S[:2,2:];P=S[2:,2:]
    return sum(H[a,b]*H[c,d]*(Q[a,c]*P[b,d]-D[a,d]*D[c,b])
               for a in range(2) for b in range(2)
               for c in range(2) for d in range(2))/M**2


def gram_from_cart(r,v):
    B=relative(r);V=relative(v)
    Y=np.concatenate((B,V),axis=1)
    return Y.T @ Y


def cart_rhs(t,z,m,G=1.0):
    r=z[:9].reshape((3,3));v=z[9:].reshape((3,3))
    a=np.zeros_like(r)
    for i,j in ((0,1),(0,2),(1,2)):
        d=r[j]-r[i]
        f=G*d/np.linalg.norm(d)**3
        a[i]+=m[j]*f
        a[j]-=m[i]*f
    return np.r_[v.ravel(),a.ravel()]


def step(S,m,h,G=1.0,rtol=1e-13,maxit=25):
    # Implicit-midpoint time-centering for w; every iterate is congruence.
    w=getw(S,G)
    for it in range(1,maxit+1):
        K=Kmat(m,w)
        Th=expm(0.5*h*K)
        Smid=Th.T@S@Th
        wnew=getw(Smid,G)
        if max(abs(wnew-w)) <= rtol*max(1.0,max(abs(wnew))):
            w=wnew
            break
        w=wnew
    else:
        raise RuntimeError('w orta-nokta sabit nokta iterasyonu yakinsamadi')
    T=expm(h*Kmat(m,w))
    Sn=T.T@S@T
    Sn=(Sn+Sn.T)/2
    return Sn,it


def integrate(S0,m,t,steps,G=1.0):
    S=S0.copy(); its=[]
    for k in range(steps):
        S,it=step(S,m,t/steps,G)
        its.append(it)
    return S,max(its),np.mean(its)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--T',type=float,default=1.0)
    args=parser.parse_args()
    m=np.array([1.,2.5,.7]);G=1.0
    r0=np.array([[-1.0,.2,.3],[1.,0.,.7],[0.,1.5,-.4]])
    v0=np.array([[.1,.4,.2],[-.2,.1,-.1],[.3,-.2,.15]])
    S0=gram_from_cart(r0,v0)
    w=getw(S0,G); A=Amat(m,w); H=Hmat(m)
    algebra=float(np.max(np.abs(A@H-H@A.T)))
    print(f'CEBIR  max |AH-H A^T| = {algebra:.3e}')
    assert algebra < 1e-12
    E0=kinetic_potential(S0,m,G);L20=L2_from_gram(S0,m)
    print(f'BASLANGIC E={E0:.15g}, L2={L20:.15g}, detS={np.linalg.det(S0):.3e}')
    z0=np.r_[r0.ravel(),v0.ravel()]
    truth=solve_ivp(cart_rhs,[0,args.T],z0,args=(m,G),rtol=1e-13,atol=1e-13,method='DOP853')
    if not truth.success:raise RuntimeError(truth.message)
    z=truth.y[:,-1];ref=gram_from_cart(z[:9].reshape(3,3),z[9:].reshape(3,3))
    prev=None
    for n in (50,100,200,400):
        S,maxit,meanit=integrate(S0,m,args.T,n,G)
        error=float(np.max(np.abs(S-ref)))
        Eerr=abs(kinetic_potential(S,m,G)-E0)
        Lerr=abs(L2_from_gram(S,m)-L20)
        determinant=abs(np.linalg.det(S))
        eigmin=np.linalg.eigvalsh(S)[0]
        ratio=(prev/error if prev else np.nan)
        print(f'step={args.T/n:.5g} n={n:4} errS={error:.4e} ratio={ratio:.3g} '
              f'dE={Eerr:.3e} dL2={Lerr:.3e} detS={determinant:.3e} eigmin={eigmin:.3e} '
              f'it(max/avg)={maxit}/{meanit:.2f}')
        assert Lerr<1e-9
        assert determinant<1e-9
        prev=error
    print('SONUC: yapisal korunum; hata orani asimptotik olarak ~4 (2. mertebe)')


if __name__=='__main__':
    main()