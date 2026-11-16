#!/usr/bin/env python3
"""
uc_cisim_tur7.py -- Tur 7: aralik (top) aritmetigiyle sertifika.  Gerekli: python-flint (Arb).   m = (1,1,1), E = 0, L = 0.

  python3 uc_cisim_tur7.py label 57 [N_tohum] [bit]     tek aci: dogrulanmis tohum + kapsamali akis + kesin etiket
  python3 uc_cisim_tur7.py check                         cizge alanlarinin kayan noktali tur-6 alanlariyla karsilastirmasi

TEOREM A (bilgisayar destekli).  Gamma(57 derece) ve Gamma(58 derece) W^u(L-) noktalari icin kacan cisim sirasiyla 1 ve 2'dir.

Kesinlik zinciri:
 (A) Isin.  Sistem 0: Z = (s, y, u_k = 1/d_k), 11 degisken, polinom alan F.  Lam t Z'(t) = F(Z(t)), Z(0) = L-, Z'(0) = V(psi).
     Z = Z^N + W;  kuyruk W, T(W)_n = (n Lam - A)^-1 [G(Z^N + W)]_n (n > N) donusumunun sabit noktasi, A = DF(L-), G = F - A.
     Norm ||x|| = max_i sum_n |x_{n,i}| t0^n.  K = 1/((N+1) Lam - ||A||_inf),  Y = ||Pi_{>N} G(Z^N)||,  g = G'nin majorani.
     K g'(xi_N + rho) <= 1/2 ve rho = 2 K Y  =>  |Z_i(t) - Z^N_i(t)| <= rho, |t| <= t0.
 (B) Akis.  Sistem 1: (0,2|1) haritasinda Levi-Civita degiskenleri (zeta, s2, omega, y2, uA, uB), 10 degisken, polinom alan
     (ikili yakin gecislerde duzgun).  Her adim: Taylor katsayilari top aritmetigiyle; kuyruk Cauchy siniriyla:
     karmasik cok-diskte |F_i| <= M_i  =>  cozum |sigma| < rho_P = min_i R/M_i icin analitik ve |a_{i,n}| <= R / rho_P^n.
 (C) Etiket.  H = 0, |s| = 1  =>  U(s) >= nu^2/2.  En kucuk iki mesafenin esitlendigi duvarlarda U <= U_E = 5/sqrt 2.
     nu azalmaz; nu > nu_E = sqrt(2 U_E) bir kez saglaninca en yakin cift sonsuza dek ayni kalir; kacan = ucuncu cisim.
"""
import sys, time, math, operator
from flint import arb, acb, arb_mat, arb_poly, ctx

# ------------------------------------------------------------------ ifade cizgesi
class Node:
    __slots__ = ("op", "a", "b", "k", "c", "idx")

    def __init__(self, op, a=None, b=None, k=None, idx=None):
        self.op, self.a, self.b, self.k, self.idx, self.c = op, a, b, k, idx, []

    def __add__(self, o): return Node("add", self, o)
    def __sub__(self, o): return Node("sub", self, o)
    def __neg__(self): return Node("scale", self, k=arb(-1))

    def __mul__(self, o):
        return Node("mul", self, o) if isinstance(o, Node) else Node("scale", self, k=arb(o))

    __rmul__ = __mul__

    def get(self, n):
        """Taylor ozyinelemesi: n. katsayi (daha dusukler hesaplanmis varsayilir)."""
        c = self.c
        if len(c) > n:
            return c[n]
        op = self.op
        if op == "mul":
            self.a.get(n); self.b.get(n)
            v = sum(map(operator.mul, self.a.c[:n + 1], reversed(self.b.c[:n + 1])), arb(0))
        elif op == "add":
            v = self.a.get(n) + self.b.get(n)
        elif op == "sub":
            v = self.a.get(n) - self.b.get(n)
        elif op == "scale":
            v = self.k * self.a.get(n)
        else:
            raise IndexError("degisken katsayisi yok")
        c.append(v)
        return v


def nodes_of(outs):
    seen, order = set(), []

    def visit(x):
        if id(x) in seen:
            return
        seen.add(id(x))
        if x.a is not None: visit(x.a)
        if x.b is not None: visit(x.b)
        order.append(x)

    for o in outs:
        visit(o)
    return order


def evaluate(order, outs, vals, mode="exact", trunc=None):
    """mode = exact: arb ya da arb_poly (trunc ile kesme);  mode = major: mutlak katsayili majorant (cikarma -> toplama)."""
    memo = {}
    for x in order:
        if x.op == "var":
            v = vals[x.idx]
        elif x.op == "mul":
            v = memo[id(x.a)] * memo[id(x.b)]
            if trunc:
                v = v.truncate(trunc)
        elif x.op == "add" or (x.op == "sub" and mode == "major"):
            v = memo[id(x.a)] + memo[id(x.b)]
        elif x.op == "sub":
            v = memo[id(x.a)] - memo[id(x.b)]
        else:
            v = memo[id(x.a)] * (abs(x.k) if mode == "major" else x.k)
        memo[id(x)] = v
    return [memo[id(o)] for o in outs]


class System:
    def __init__(self, nvar, build):
        self.n = nvar
        self.vars = [Node("var", idx=i) for i in range(nvar)]
        self.outs = build(self.vars)
        self.order = nodes_of(self.outs)

    def F(self, Z, trunc=None):
        return evaluate(self.order, self.outs, Z, "exact", trunc)

    def major(self, Zabs):
        return evaluate(self.order, self.outs, [arb_poly([z, arb(1)]) for z in Zabs], "major")

    def taylor(self, Z0, N):
        for x in self.order:
            x.c = []
        for i, v in enumerate(self.vars):
            v.c = [Z0[i]]
        for n in range(N):
            f = [o.get(n) for o in self.outs]
            for i, v in enumerate(self.vars):
                v.c.append(f[i] / (n + 1))
        return [v.c for v in self.vars]


def up(x):
    return x.abs_upper()


# ------------------------------------------------------------------ sistem 0: (s, y, u), Jacobi (0,1|2), m = 1
def build0(Z):
    r2 = arb(2).sqrt()
    AL, BE = [r2, r2 / 2, -r2 / 2], [arb(0), (arb(3) / 2).sqrt(), (arb(3) / 2).sqrt()]
    s, y, u = Z[:4], Z[4:8], Z[8:]
    nu = s[0] * y[0] + s[1] * y[1] + s[2] * y[2] + s[3] * y[3]
    ds = [y[a] - nu * s[a] for a in range(4)]
    dy = [(nu * y[a]) * (arb(1) / 2) for a in range(4)]
    du = []
    for k in range(3):
        D = [s[0] * AL[k] + s[2] * BE[k], s[1] * AL[k] + s[3] * BE[k]]
        Dp = [ds[0] * AL[k] + ds[2] * BE[k], ds[1] * AL[k] + ds[3] * BE[k]]
        u3 = (u[k] * u[k]) * u[k]
        w = [u3 * D[0], u3 * D[1]]
        dy[0], dy[1] = dy[0] - w[0] * AL[k], dy[1] - w[1] * AL[k]
        dy[2], dy[3] = dy[2] - w[0] * BE[k], dy[3] - w[1] * BE[k]
        du.append(-(u3 * (D[0] * Dp[0] + D[1] * Dp[1])))
    return ds + dy + du


# ------------------------------------------------------------------ sistem 1: Levi-Civita, harita (0,2|1), m = 1, E = 0
def lc_consts():
    return dict(c2=(arb(3) / 2).sqrt(), a=arb(2).sqrt() / 2, r1=(arb(1) / 2).sqrt(), r2=(arb(2) / 3).sqrt())


def build1(Z):
    K = lc_consts()
    c2, a = K["c2"], K["a"]
    zr, zi, sr, si, wr, wi, yr, yi, uA, uB = Z
    zz_r, zz_i = zr * zr, zi * zi
    z2 = zz_r + zz_i
    nu = zr * wr + zi * wi + sr * yr + si * yi
    P = nu * z2
    dzr, dzi = (wr - P * zr) * (arb(1) / 2), (wi - P * zi) * (arb(1) / 2)
    s1r, s1i = zz_r - zz_i, (zr * zi) * 2
    Ar, Ai = sr * c2 + s1r * a, si * c2 + s1i * a
    Br, Bi = sr * c2 - s1r * a, si * c2 - s1i * a
    uA3, uB3 = (uA * uA) * uA, (uB * uB) * uB
    fAr, fAi, fBr, fBi = Ar * uA3, Ai * uA3, Br * uB3, Bi * uB3
    g1r, g1i = (fBr - fAr) * a, (fBi - fAi) * a
    g2r, g2i = -((fAr + fBr) * c2), -((fAi + fBi) * c2)
    W = uA + uB - (yr * yr + yi * yi) * (arb(1) / 2)
    dwr = zr * W + z2 * (zr * g1r + zi * g1i)
    dwi = zi * W + z2 * (zr * g1i - zi * g1r)
    dsr, dsi = z2 * yr - P * sr, z2 * yi - P * si
    dyr, dyi = z2 * g2r + (P * yr) * (arb(1) / 2), z2 * g2i + (P * yi) * (arb(1) / 2)
    qr, qi = zr * dzr - zi * dzi, zr * dzi + zi * dzr                 # zeta zeta'
    dAr, dAi = dsr * c2 + qr * (2 * a), dsi * c2 + qi * (2 * a)
    dBr, dBi = dsr * c2 - qr * (2 * a), dsi * c2 - qi * (2 * a)
    duA = -(uA3 * (Ar * dAr + Ai * dAi))
    duB = -(uB3 * (Br * dBr + Bi * dBi))
    return [dzr, dzi, dsr, dsi, dwr, dwi, dyr, dyi, duA, duB]


def to_lc(Z0):
    """Sistem 0 topu (s, y) -> sistem 1 topu. Jacobi (0,1|2) -> (0,2|1) dondurmesi + zeta = sqrt(s1)."""
    K = lc_consts()
    r1, r2, c2, a = K["r1"], K["r2"], K["c2"], K["a"]
    rt2 = arb(2).sqrt()

    def rot(q):
        x = [rt2 * q[0], rt2 * q[1]]
        X = [c2 * q[2], c2 * q[3]]
        return [r1 * (X[0] + x[0] / 2), r1 * (X[1] + x[1] / 2), r2 * (arb(3) / 4 * x[0] - X[0] / 2), r2 * (arb(3) / 4 * x[1] - X[1] / 2)]

    s, y = rot(Z0[:4]), rot(Z0[4:8])
    ze = acb(s[0], s[1]).sqrt()
    om = ze.conjugate() * acb(y[0], y[1]) if hasattr(ze, "conjugate") else acb(ze.real, -ze.imag) * acb(y[0], y[1])
    s1 = ze * ze
    A = [c2 * s[2] + a * s1.real, c2 * s[3] + a * s1.imag]
    B = [c2 * s[2] - a * s1.real, c2 * s[3] - a * s1.imag]
    uA, uB = 1 / (A[0] * A[0] + A[1] * A[1]).sqrt(), 1 / (B[0] * B[0] + B[1] * B[1]).sqrt()
    return [ze.real, ze.imag, s[2], s[3], om.real, om.imag, y[2], y[3], uA, uB]


def lc_observe(Z):
    K = lc_consts()
    c2, a = K["c2"], K["a"]
    zr, zi, sr, si, wr, wi, yr, yi, uA, uB = Z
    z2 = zr * zr + zi * zi
    nu = zr * wr + zi * wi + sr * yr + si * yi
    s1r, s1i = zr * zr - zi * zi, 2 * zr * zi
    A2 = (c2 * sr + a * s1r) ** 2 + (c2 * si + a * s1i) ** 2
    B2 = (c2 * sr - a * s1r) ** 2 + (c2 * si - a * s1i) ** 2
    d2 = {2: 2 * z2 * z2, 3: A2, 1: B2}                # anahtar = cift disinda kalan cisim: (1,3)->2, (1,2)->3, (2,3)->1
    return nu, d2, z2


# ------------------------------------------------------------------ (A) dogrulanmis isin
def lagrange():
    r3 = arb(3).sqrt()
    pos = [(arb(0), 1 / r3), (-arb(1) / 2, -1 / (2 * r3)), (arb(1) / 2, -1 / (2 * r3))]
    r1, r2 = (arb(1) / 2).sqrt(), (arb(2) / 3).sqrt()

    def to_q(p):
        c = [(p[0][i] + p[1][i]) / 2 for i in range(2)]
        return [r1 * (p[1][0] - p[0][0]), r1 * (p[1][1] - p[0][1]), r2 * (p[2][0] - c[0]), r2 * (p[2][1] - c[1])]

    sL = to_q(pos)
    Js = [-sL[1], sL[0], -sL[3], sL[2]]
    dot = lambda u, v: sum((u[i] * v[i] for i in range(4)), arb(0))
    base = []
    for d in ((arb(0), arb(1)), (arb(-1), arb(0))):           # cisim 1: disari / tegetsel
        q = to_q([d, (arb(0), arb(0)), (arb(0), arb(0))])
        al, be = dot(q, sL), dot(q, Js)
        q = [q[i] - al * sL[i] - be * Js[i] for i in range(4)]
        n = dot(q, q).sqrt()
        base.append([x / n for x in q])
    vL = arb(6).sqrt()
    return sL + [-vL * x for x in sL] + [arb(1)] * 3, sL, base, vL


def validated_seed(S0, psi_deg, t0, N, verbose=True):
    Zeq, sL, base, vL = lagrange()
    Lam = vL / 4 * (1 + arb(13).sqrt())
    ps = arb(psi_deg) * arb.pi() / 180
    sig = [ps.cos() * base[0][i] + ps.sin() * base[1][i] for i in range(4)]
    r2 = arb(2).sqrt()
    AL, BE = [r2, r2 / 2, -r2 / 2], [arb(0), (arb(3) / 2).sqrt(), (arb(3) / 2).sqrt()]
    du = []
    for k in range(3):
        D = [AL[k] * sL[0] + BE[k] * sL[2], AL[k] * sL[1] + BE[k] * sL[3]]
        Ds = [AL[k] * sig[0] + BE[k] * sig[2], AL[k] * sig[1] + BE[k] * sig[3]]
        du.append(-(D[0] * Ds[0] + D[1] * Ds[1]))
    V = sig + [(Lam - vL) * x for x in sig] + du
    cols = []
    for k in range(11):
        f = S0.F([arb_poly([Zeq[i], arb(1 if i == k else 0)]) for i in range(11)], trunc=2)
        cols.append([(x[1] if x.length() > 1 else arb(0)) for x in f])
    A = [[cols[j][i] for j in range(11)] for i in range(11)]
    eq_res = max(up(x) for x in S0.F(Zeq))
    ev_res = max(up(sum((A[i][j] * V[j] for j in range(11)), arb(0)) - Lam * V[i]) for i in range(11))
    coef = [[Zeq[i], V[i]] for i in range(11)]
    for n in range(2, N + 1):
        Fn = S0.F([arb_poly(coef[i]) for i in range(11)], trunc=n + 1)
        rhs = arb_mat(11, 1, [(Fn[i][n] if Fn[i].length() > n else arb(0)) for i in range(11)])
        Mat = arb_mat(11, 11, [(n * Lam if i == j else arb(0)) - A[i][j] for i in range(11) for j in range(11)])
        sol = Mat.solve(rhs)
        for i in range(11):
            coef[i].append(sol[i, 0])
    t0a = arb(t0)
    normA = max(sum((up(A[i][j]) for j in range(11)), arb(0)) for i in range(11))
    den = (N + 1) * Lam - normA
    assert den > 0
    K = up(1 / den)
    xiN = max(sum((up(coef[i][n]) * t0a ** n for n in range(1, N + 1)), arb(0)) for i in range(11))
    Q = S0.major([up(z) for z in Zeq])
    Ff = S0.F([arb_poly(coef[i]) for i in range(11)])
    Y = arb(0)
    for i in range(11):
        cf = Ff[i].coeffs()
        Y = Y.max(sum((up(cf[n]) * t0a ** n for n in range(N + 1, len(cf))), arb(0)))
    Y = up(Y)
    rho = up(2 * K * Y)
    xi = up(xiN + rho)
    gp = arb(0)
    for i in range(11):
        cf = Q[i].coeffs()
        gp = gp.max(sum((k * up(cf[k]) * xi ** (k - 1) for k in range(2, len(cf))), arb(0)))
    kap = up(K * gp)
    ok = bool(kap < arb(1) / 2)
    Z = [sum((coef[i][n] * t0a ** n for n in range(N + 1)), arb(0)) + arb(0, rho) for i in range(11)]
    info = dict(ok=ok, contraction=float(kap), rho=float(rho), Y=float(Y), K=float(K), normA=float(normA), xiN=float(xiN),
                eq_res=float(eq_res), ev_res=float(ev_res))
    if verbose:
        print(f"   [A] psi = {psi_deg}, t0 = {t0}, N = {N}, {ctx.prec} bit: ||A||_inf <= {info['normA']:.3f}, K <= {info['K']:.3e},"
              f" xi_N <= {info['xiN']:.4f}, Y <= {info['Y']:.3e}", flush=True)
        print(f"       daralma K g'(xi_N + rho) <= {info['contraction']:.4f} (< 1/2: {ok});  kuyruk siniri rho = {info['rho']:.3e};"
              f"  |F(L-)| <= {info['eq_res']:.0e}, |(A - Lam) V| <= {info['ev_res']:.0e}", flush=True)
    return Z, info


# ------------------------------------------------------------------ (B) kapsamali akis
def picard(S, Z, Rs=(0.005, 0.01, 0.02, 0.04, 0.08, 0.16)):
    F0 = [up(f) for f in S.F(Z)]
    Q = [q.coeffs() for q in S.major([up(z) for z in Z])]
    best = None
    for R in Rs:
        Ra, rho = arb(R), None
        for i in range(S.n):
            M = F0[i] + sum((up(Q[i][k]) * Ra ** k for k in range(1, len(Q[i]))), arb(0))
            r = (Ra / M).lower()
            rho = r if rho is None else rho.min(r)
        if best is None or float(rho) > float(best[1]):
            best = (R, rho.lower())
    return best


def horner(c, h):
    v = arb(0)
    for x in reversed(c):
        v = v * h + x
    return v


def flow_to_decision(S1, Z, theta=0.2, margin=1e-3, max_steps=40000, log=None):
    nuE = (5 * arb(2).sqrt()).sqrt()
    target = nuE + arb(margin)
    sig, steps, Nmax, t_start = 0.0, 0, 0, time.time()
    rad0 = max(float(z.rad()) for z in Z)
    hmin = math.inf
    while steps < max_steps:
        nu, d2, z2 = lc_observe(Z)
        if nu > target:
            break
        rad = max(float(z.rad()) for z in Z)
        if not rad < 1e-4:
            return dict(ok=False, why="top cok genisledi", sigma=sig, steps=steps, rad=rad, nu=float(nu.mid()))
        R, rhoP = picard(S1, Z)
        h = float(rhoP) * theta * (1 - 1e-9)
        ha = arb(h)
        th = up(ha / rhoP)
        N = max(8, int(math.ceil(math.log(R / (0.02 * (1 - theta) * max(rad, 1e-300))) / math.log(1 / theta))))
        co = S1.taylor(Z, N)
        tail = up(arb(R) * th ** (N + 1) / (1 - th))
        Z = [horner(co[i], ha) + arb(0, tail) for i in range(S1.n)]
        sig, steps, Nmax, hmin = sig + h, steps + 1, max(Nmax, N), min(hmin, h)
        if log and steps % log == 0:
            print(f"       adim {steps}: sigma = {sig:.4f}, nu = {float(nu.mid()):+.5f}, |zeta|^2 = {float(z2.mid()):.4f}, R = {R},"
                  f" rho_P = {float(rhoP):.5f}, N = {N}, yaricap = {max(float(z.rad()) for z in Z):.2e}, sure = {time.time() - t_start:.0f} s", flush=True)
    else:
        return dict(ok=False, why="adim siniri", sigma=sig, steps=steps)
    close = [k for k in d2 if all(d2[k] < d2[j] for j in d2 if j != k)]
    zr, zi, sr, si, wr, wi, yr, yi, uA, uB = Z
    Gam = (wr * wr + wi * wi) / 2 - (arb(1) / 2).sqrt() - z2 * (uA + uB - (yr * yr + yi * yi) / 2)
    return dict(ok=len(close) == 1, label=close[0] if close else None, sigma=sig, steps=steps, Nmax=Nmax, hmin=hmin, rad0=rad0,
                rad=max(float(z.rad()) for z in Z), nu_lower=float(nu.lower()), nuE=float(nuE.mid()),
                dist={k: float(v.sqrt().mid()) for k, v in d2.items()}, Gamma=float(up(Gam)),
                norm=float(up(z2 * z2 + sr * sr + si * si - 1)), time=time.time() - t_start, Zmid=[float(z.mid()) for z in Z])


def certify(psi, Nseed=150, bits=640, t0=0.1, log=200):
    ctx.prec = bits
    S0, S1 = System(11, build0), System(10, build1)
    Z, info = validated_seed(S0, psi, t0, Nseed)
    if not info["ok"]:
        return dict(ok=False, why="tohum daralmasi saglanmadi")
    out = flow_to_decision(S1, to_lc(Z), log=log)
    out["seed"] = info
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"
    if mode == "label":
        psi = float(sys.argv[2])
        out = certify(psi, int(sys.argv[3]) if len(sys.argv) > 3 else 150, int(sys.argv[4]) if len(sys.argv) > 4 else 640)
        print("   [B] " + ", ".join(f"{k} = {v}" for k, v in out.items() if k not in ("seed", "Zmid", "dist")))
        if out.get("ok"):
            print(f"   [C] psi = {psi}: nu >= {out['nu_lower']:.6f} > nu_E = {out['nuE']:.6f};  mesafeler (kacan adayina gore) = {out['dist']}")
            print(f"   => KESIN: Gamma({psi} derece) icin kacan cisim = {out['label']}")
    elif mode == "check":
        import numpy as np
        sys.path.insert(0, ".")
        import uc_cisim_tur6 as T6
        ctx.prec = 200
        S0, S1 = System(11, build0), System(10, build1)
        W = T6.Manifold((1, 1, 1), 24)
        a = (0.1 * math.cos(math.radians(57)), 0.1 * math.sin(math.radians(57)))
        Zf = W(a)
        F0 = [float(x.mid()) for x in S0.F([arb(float(z)) for z in Zf])]
        print("   sistem 0: max |F_cizge - F_tur6| =", np.abs(np.array(F0) - W.P.Ffloat(Zf)).max())
        Z, info = validated_seed(S0, 57, 0.1, 60)
        print("   tohum: max |orta - tur6 P(a)| =", max(abs(float(Z[i].mid()) - Zf[i]) for i in range(11)))
        lc = T6.LC(np.ones(3), 1)
        Zl = to_lc([arb(float(z)) for z in Zf])
        r, v = T6.Jac(np.ones(3)).to_r(Zf[:4]), T6.Jac(np.ones(3)).to_r(Zf[4:8])
        X = lc.from_cart(r, v, 0.0, 0.0)
        sgn = 1.0 if abs(float(Zl[0].mid()) - X[0]) < 1e-9 else -1.0
        print("   harita donusumu: max |to_lc - tur6| =", max(abs(float(Zl[i].mid()) - sgn * X[i]) for i in (0, 1, 4, 5)),
              max(abs(float(Zl[i].mid()) - X[i]) for i in (2, 3, 6, 7)))
        F1 = [float(x.mid()) for x in S1.F(Zl)]
        Xs = [sgn * X[0], sgn * X[1], X[2], X[3], sgn * X[4], sgn * X[5], X[6], X[7], 0.0, 0.0]
        ref = lc.rhs(0.0, Xs)
        print("   sistem 1: max |F_cizge - LC.rhs(tur6)| =", max(abs(F1[i] - ref[i]) for i in range(8)))
        nu, d2, z2 = lc_observe(Zl)
        rr = T6.Jac(np.ones(3)).to_r(Zf[:4])
        dd = {2: np.linalg.norm(rr[2] - rr[0]), 3: np.linalg.norm(rr[1] - rr[0]), 1: np.linalg.norm(rr[2] - rr[1])}
        print("   gozlemler: nu farki =", abs(float(nu.mid()) - Zf[:4] @ Zf[4:8]), " mesafe farki =", max(abs(float(d2[k].sqrt().mid()) - dd[k]) for k in dd))
        co = S1.taylor(Zl, 12)
        h = 1e-3
        Zh = [float(horner(co[i], arb(h)).mid()) for i in range(10)]
        from scipy.integrate import solve_ivp
        sol = solve_ivp(lambda t, x: [float(q.mid()) for q in S1.F([arb(float(c)) for c in x])], [0, h], [float(z.mid()) for z in Zl], rtol=1e-13, atol=1e-15, method="DOP853")
        print("   Taylor adimi vs DOP853 (h = 1e-3): max fark =", max(abs(Zh[i] - sol.y[i, -1]) for i in range(10)))