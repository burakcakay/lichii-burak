"""G03/G04 independent numerical and symbolic smoke tests.

Passing this program does not establish global integrability or novelty.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

from g03_quartic import (
    mass_coefficients, encode, unpack, force_matrix,
    quartic_rhs, newton_rhs, angular_momentum_sq,
    oriented_minors, taylor_coefficients, taylor_evaluate,
)

rng = np.random.default_rng(20261121)
max_gram = max_h = max_energy = max_j = max_alg = max_taylor = 0.0
for case in range(5):
    masses = rng.uniform(0.6, 2.0, 3)
    c, strengths = mass_coefficients(masses)
    v = rng.normal(size=(4, 3)) * 0.35
    v[0] = [0.6, 0.15, -0.2]
    v[1] = [0.25, 1.1, 0.2]
    if case == 4:
        v[:, 2] = 0
        assert np.max(abs(oriented_minors(v))) < 1e-14

    z0 = encode(v, c)
    energy0 = 0.5 * np.sum(v[2:]**2) - np.dot(strengths, z0[10:13])
    mom0 = np.cross(v[0], v[2]) + np.cross(v[1], v[3])
    mom2_0 = mom0 @ mom0

    invariant = solve_ivp(quartic_rhs, (0, 0.12), z0, args=(c, strengths),
                          method="DOP853", rtol=1e-12, atol=1e-12,
                          dense_output=True)
    reference = solve_ivp(newton_rhs, (0, 0.12), v.ravel(), args=(c, strengths),
                          method="DOP853", rtol=1e-12, atol=1e-12,
                          dense_output=True)
    assert invariant.success and reference.success
    for t in np.linspace(0.0, 0.12, 20):
        z = invariant.sol(t)
        vv = reference.sol(t).reshape(4, 3)
        q, mixed, p, h = unpack(z)
        gram = np.block([[q, mixed], [mixed.T, p]])
        max_gram = max(max_gram, np.max(abs(gram - vv @ vv.T)))
        max_h = max(max_h,
                    np.max(abs(h - 1.0 / np.linalg.norm(c @ vv[:2], axis=1))))
        max_alg = max(max_alg, np.max(abs(
            h*h * np.einsum("ia,ab,ib->i", c, q, c) - 1)))
        max_energy = max(max_energy,
                         abs(0.5*np.trace(p) - np.dot(strengths,h) - energy0))
        max_j = max(max_j, abs(angular_momentum_sq(q,mixed,p) - mom2_0))

    zt = taylor_evaluate(taylor_coefficients(z0, c, strengths, 12), 0.03)
    max_taylor = max(max_taylor, np.max(abs(zt - invariant.sol(0.03))))

assert max_gram < 3e-10, max_gram
assert max_h < 3e-10, max_h
assert max_energy < 1e-10, max_energy
assert max_j < 3e-10, max_j
assert max_alg < 3e-10, max_alg
assert max_taylor < 3e-10, max_taylor

O=sp.zeros(2)
I=sp.eye(2)
def block(a,b,c,d): return a.row_join(b).col_join(c.row_join(d))
def L(k): return block(O,O,k,O)
def comm(a,b): return a*b-b*a
sym_basis=[
    sp.Matrix([[1,0],[0,0]]),
    sp.Matrix([[0,0],[0,1]]),
    sp.Matrix([[0,1],[1,0]]),
]
basis=[block(O,I,O,O)]+[L(x) for x in sym_basis]
def rank(mats):
    return sp.Matrix.hstack(*[sp.Matrix(x).reshape(16,1) for x in mats]).rank()
for _ in range(3):
    old=list(basis)
    for a in old:
        for b in old:
            ab=comm(a,b)
            if ab!=sp.zeros(4) and rank(basis+[ab])>len(basis):
                basis.append(ab)
    if len(basis)==10:
        break
assert len(basis)==10
symplectic=block(O,I,-I,O)
assert all(a.T*symplectic+symplectic*a==sp.zeros(4) for a in basis)

print("PASS G03: max Gram", f"{max_gram:.4e}",
      "max h", f"{max_h:.4e}", "max energy",f"{max_energy:.4e}",
      "max J²", f"{max_j:.4e}", "constraints",f"{max_alg:.4e}")
print("PASS G03: degree-12 Taylor", f"{max_taylor:.4e}")
print("PASS G04: full sp(4,R) Lie rank", len(basis))
