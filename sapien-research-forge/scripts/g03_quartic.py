"""G03: quartic physical-time invariant flow for Newtonian three bodies.

Research lead: Burak Cakay (@burakcakay). This implements an exact
representation (not a general closed-form solution).
"""
from __future__ import annotations

import math
import numpy as np


def mass_coefficients(masses):
    m1, m2, m3 = map(float, masses)
    if min(m1, m2, m3) <= 0:
        raise ValueError("all masses must be positive")
    m12 = m1 + m2
    mu1 = m1 * m2 / m12
    mu2 = m3 * m12 / (m12 + m3)
    c = np.array([
        [1.0 / math.sqrt(mu1), 0.0],
        [m2 / (m12 * math.sqrt(mu1)), 1.0 / math.sqrt(mu2)],
        [-m1 / (m12 * math.sqrt(mu1)), 1.0 / math.sqrt(mu2)],
    ])
    return c, np.array([m1 * m2, m1 * m3, m2 * m3])


def encode(v, c):
    """Convert four rows (q1,q2,p1,p2) of three-vectors to 13 invariants."""
    gram = v @ v.T
    distances = np.linalg.norm(c @ v[:2], axis=1)
    if np.any(distances <= 0):
        raise ValueError("binary collision")
    return np.r_[
        gram[0, 0], gram[0, 1], gram[1, 1],
        gram[:2, 2:].ravel(),
        gram[2, 2], gram[2, 3], gram[3, 3],
        1.0 / distances,
    ]


def unpack(z):
    q = np.array([[z[0], z[1]], [z[1], z[2]]])
    mixed = z[3:7].reshape(2, 2)
    p = np.array([[z[7], z[8]], [z[8], z[9]]])
    h = z[10:13]
    return q, mixed, p, h


def force_matrix(h, c, strengths):
    return np.einsum("i,ia,ib->ab", strengths * h**3, c, c)


def quartic_rhs(t, z, c, strengths):
    """Polynomial in z of total degree <=4; physical Newtonian t."""
    q, mixed, p, h = unpack(z)
    k = force_matrix(h, c, strengths)
    dq = mixed + mixed.T
    dm = p - q @ k
    dp = -k @ mixed - mixed.T @ k
    dh = -(h**3) * np.einsum("ia,ab,ib->i", c, mixed, c)
    return np.r_[
        dq[0, 0], dq[0, 1], dq[1, 1],
        dm.ravel(), dp[0, 0], dp[0, 1], dp[1, 1], dh,
    ]


def newton_rhs(t, vector_state, c, strengths):
    """Independent 12D reference, still in COM Jacobi coordinates."""
    v = vector_state.reshape(4, 3)
    d = np.linalg.norm(c @ v[:2], axis=1)
    k = force_matrix(1.0 / d, c, strengths)
    return np.vstack([v[2:], -k @ v[:2]]).ravel()


def angular_momentum_sq(q, mixed, p):
    return (
        q[0, 0] * p[0, 0] - mixed[0, 0] ** 2
        + q[1, 1] * p[1, 1] - mixed[1, 1] ** 2
        + 2 * (q[0, 1] * p[0, 1] - mixed[0, 1] * mixed[1, 0])
    )


def oriented_minors(v):
    return np.array([
        (-1)**i * np.linalg.det(np.delete(v, i, axis=0))
        for i in range(4)
    ])


def taylor_coefficients(z0, c, strengths, order):
    """Exact algebraic coefficient recurrence computed in floating point."""
    q = np.zeros((order + 1, 2, 2))
    mixed = np.zeros_like(q)
    p = np.zeros_like(q)
    h = np.zeros((order + 1, 3))
    k = np.zeros_like(q)
    q[0], mixed[0], p[0], h[0] = unpack(z0)

    for n in range(order):
        h2 = np.zeros((n + 1, 3))
        h3 = np.zeros((n + 1, 3))
        for j in range(n + 1):
            h2[j] = sum((h[i] * h[j-i] for i in range(j + 1)), np.zeros(3))
            h3[j] = sum((h2[i] * h[j-i] for i in range(j + 1)), np.zeros(3))
        k[n] = np.einsum("i,ia,ib->ab", strengths * h3[n], c, c)
        dq = mixed[n] + mixed[n].T
        dm = p[n] - sum(
            (q[j] @ k[n-j] for j in range(n + 1)), np.zeros((2, 2)))
        dp = -sum(
            (k[j] @ mixed[n-j] + mixed[n-j].T @ k[j]
             for j in range(n + 1)), np.zeros((2, 2)))
        dh = -sum(
            (h3[j] * np.einsum("ia,ab,ib->i", c, mixed[n-j], c)
             for j in range(n + 1)), np.zeros(3))
        q[n+1], mixed[n+1], p[n+1], h[n+1] = (
            dq/(n+1), dm/(n+1), dp/(n+1), dh/(n+1))
    return list(zip(q, mixed, p, h))


def taylor_evaluate(coefficients, t):
    q = np.zeros((2, 2)); mixed = np.zeros((2, 2))
    p = np.zeros((2, 2)); h = np.zeros(3)
    for cq, cm, cp, ch in reversed(coefficients):
        q = q * t + cq; mixed = mixed * t + cm
        p = p * t + cp; h = h * t + ch
    return np.r_[q[0,0], q[0,1], q[1,1], mixed.ravel(),
                 p[0,0], p[0,1], p[1,1], h]
