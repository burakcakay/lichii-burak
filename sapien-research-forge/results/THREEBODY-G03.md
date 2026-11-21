# G03 — Universal quartic Gram flow for the Newtonian three-body system

**Status: ANALYTIC_PROOF_CANDIDATE.** Re-derived and independently numerically cross-checked, but not independently peer reviewed; **no literature-priority claim**. Polynomializing inverse distances is known in N-body research (see related work at end). The specific combined Gram/chirality/quartic/self-consistent symplectic formulation is a proposed research architecture whose novelty must be checked.

## A. Domain and mass coefficients

Take arbitrary m1,m2,m3>0, G_Newton=1, physical COM frame, spatial dimension 3; planar states are included. Let m12=m1+m2, Mtot=m12+m3, μ1=m1m2/m12, μ2=m3m12/Mtot, and
\[
q_1=\sqrt{\mu_1}(r_2-r_1),\quad
q_2=\sqrt{\mu_2}\left(r_3-\frac{m_1r_1+m_2r_2}{m_{12}}\right),\quad
p_i=\dot q_i.
\]
Set three constant two-vectors
\[
c_1=(1/\sqrt{\mu_1},0),\
c_2=(m_2/(m_{12}\sqrt{\mu_1}),1/\sqrt{\mu_2}),\
c_3=(-m_1/(m_{12}\sqrt{\mu_1}),1/\sqrt{\mu_2}),
\]
and g=(m1m2,m1m3,m2m3). Physical pair distances are d_i=|c_{i1}q_1+c_{i2}q_2|.

## B. Exact 13-variable quartic equations

Define \(V=(q_1,q_2,p_1,p_2)^\top\in\mathbb R^{4\times3}\),
\[
\Gamma=VV^\top=\begin{pmatrix}Q&M\\M^\top&P\end{pmatrix},
\]
with Q,P symmetric 2×2, M arbitrary 2×2, and auxiliary inverse distances \(h_i=1/d_i>0\).
Let
\[
K(h)=\sum_{i=1}^3 g_i h_i^3 c_i c_i^\top.
\]
The **physical-time**, degree-at-most-four polynomial flow is
\[
\boxed{
\dot Q=M+M^\top,\quad
\dot M=P-QK,\quad
\dot P=-KM-M^\top K,\quad
\dot h_i=-h_i^3(c_i^\top M c_i).
}
\]
The 13 formal variables are constrained by
\[
\boxed{h_i^2(c_i^\top Qc_i)=1,\qquad i=1,2,3.}
\]
The constraints are invariant: setting \(C_i=h_i^2(c_i^\top Qc_i)\),
\[
\dot C_i=2h_i^2(c_i^\top Mc_i)(1-C_i).
\]
Thus exact initially consistent data remain consistent.

**Proof of equivalence.** The physical Newton equations in mass-weighted Jacobi coordinates are \(\dot q=p,\ \dot p=-K q\), since \(U=\sum g_i/d_i\). Differentiating \(Q=q q^\top,\ M=q p^\top,\ P=p p^\top\) gives the first three equations. Differentiating \(h_i=(c_i^\top Q c_i)^{-1/2}\) gives the fourth. Conversely, for consistent collision-free initial data, standard uniqueness plus reconstruction of a lifted V from Γ proves the Gram solution agrees with the physical Newton flow. This is an **exact change of representation**, not new global integrability.

Energy and squared total angular momentum are polynomials in lifted coordinates:
\[
E=\frac12\operatorname{tr}P-\sum_i g_i h_i,
\]
\[
\ell^2=Q_{11}P_{11}-M_{11}^2+
Q_{22}P_{22}-M_{22}^2+
2(Q_{12}P_{12}-M_{12}M_{21}).
\]
For 3-dimensional position-velocity data, rank Γ ≤3, so det Γ=0. On the regular fixed E and \(\ell>0\) stratum this yields the expected \(13-3-1-1-1=7\) quotient dimension; the seven-variable count does not make all strata smooth.

## C. Spatial handedness (SO(3), not O(3))

The Gram matrix only determines V up to O(3). Retain four oriented triple minors
\[
\chi_i=(-1)^{i-1}\det(V_{\widehat i})\quad(1\le i\le4).
\]
Set \(A=\begin{psmallmatrix}0&I\\-K&0\end{psmallmatrix}\), so \(\dot V=AV\) and \(\mathrm{tr}A=0\). By the third exterior power representation,
\[
\boxed{\dot\chi=-A^\top\chi,\quad
\Gamma\chi=0,\quad
\chi\chi^\top=\operatorname{adj}(\Gamma).}
\]
Both the exact planar condition \(\chi=0\) (rank V≤2) and the nonplanar condition \(\chi\neq0\) are dynamically invariant throughout every smooth collision-free interval.

## D. Self-consistent symplectic transport

\(A^\top\mathbb J+\mathbb J A=0\) for \(\mathbb J=\begin{psmallmatrix}0&I\\-I&0\end{psmallmatrix}\).
Therefore the fundamental transfer B satisfying \(\dot B=A(\Gamma(t))B,\ B(0)=I\) lies in Sp(4,R), and
\[
\boxed{\Gamma(t)=B(t)\Gamma(0)B(t)^\top,\quad
\chi(t)=B(t)^{-\top}\chi(0).}
\]
This is a nonlinear *self-consistent* transfer equation, not an independently computable matrix exponential.

## E. Exact Taylor recurrence

Write \(Z(t)=\sum_{n\ge0}Z_n t^n\) for the 13D polynomial state \(\dot Z=F_4(Z)\). Then
\[
\boxed{Z_{n+1}=\frac1{n+1}[t^n]F_4\left(\sum_{k=0}^nZ_kt^k\right).}
\]
Since deg F4≤4, each coefficient uses finite convolutions of previous coefficients; **no nonlinear solve** is necessary at each Taylor degree. This is a local representation only; certified complex-time analytic radii and adaptive continuation are not proved in this record.

### Measured results (Python 3.13, NumPy 2.3.5, SciPy 1.17.0)

Five deterministic noncollision examples with positive unequal masses, four spatial and one planar, direct Cartesian Jacobi DOP853 versus 13-variable quartic DOP853, t∈[0,0.12]:
- max Gram matrix discrepancy 1.1456e-11
- max inverse-distance discrepancy 4.830e-12
- max energy discrepancy 1.776e-15
- max squared angular momentum discrepancy 6.395e-12
- max reciprocal-distance constraint discrepancy 1.426e-11
- median quartic runtime ~1.77× **slower** than direct in the small Python benchmark. No performance superiority claimed.

12th degree Taylor recurrence at t=0.03: max absolute discrepancy ≤1.6e-13 across the five examples versus high-accuracy DOP853, limited by double precision. At t=0.1, one case needs 24th order to reach ~1e-14. No validated global error certificate is claimed.

## Limitations and literature

- \(h_i\to\infty\) near a physical binary collision; this is **not collision regularization**. Triple collision and collision continuation remain separate.
- Positive inverse distances, Gram rank, energy and angular momentum constrain the algebraic invariant variety.
- The reduction does not create new commuting constants or a finite closed-form global solution.
- Related: [Polynomialization of the N-body system using inverse distances](https://www.cfm.brown.edu/people/dobrush/am34/python/intro/picard.html); [On the Quadratization of the Integrals for the Many-Body Problem](https://www.mdpi.com/2227-7390/9/24/3208); [Natural dynamical reduction of the three-body problem](https://pmc.ncbi.nlm.nih.gov/articles/PMC10184114/).
