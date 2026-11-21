# T21 — Exchange-reflection symmetry cancels the quadratic resonant coefficient

**Status: ANALYTIC_PROOF_CANDIDATE.** This is an exact analytic symmetry argument with symbolic Hessian verification; independent mathematical review and literature priority check are still pending.

## The claim

Let \(m_1=m_2=1,\ m_3=M>1\) and consider the planar, zero-angular-momentum, zero-energy collapsing equilateral Lagrange rest point \(L_-\). At the algebraic mass ratio
\[
\sigma_*=(1064-32\sqrt{129})/5625,
\qquad M_*=12.38418125935020520817\ldots,
\]
the reduced shape eigenvalues obey \(\lambda_1=2\lambda_2\) and \(\lambda_4=\lambda_2+\lambda_3\). The two possible quadratic Poincaré–Dulac resonance terms on the four-dimensional \(M_0\) are
\[
c_{\rm res}\,\xi_2^2\partial_{\xi_1},
\qquad d_{\rm res}\,\xi_2\xi_3\partial_{\xi_4}.
\]

**Conclusion:** For \(M=M_*\), both coefficients vanish **exactly by symmetry**:
\[
\boxed{c_{\rm res}=d_{\rm res}=\kappa=0,\qquad\kappa=c_{\rm res}w_+.}
\]
Thus the two *quadratic* resonance channels cannot force a logarithmic term in this exchange-symmetric mass family. No claim about other nonresonant mechanisms or all higher-order effects is made.

## Explicit involution and eigendirections

Use normalized mass-weighted planar Jacobi position
\(s_0=(a,0,0,b)\in S^3\), with
\[
\mu_2=\frac{2M}{M+2},\quad
\rho_0=\sqrt{\frac{2M+1}{M+2}},\quad
a=\frac{1}{\sqrt2\,\rho_0},\quad
b=\frac{\sqrt{3\mu_2}}{2\rho_0},\quad
a^2+b^2=1.
\]
The first-two-body exchange composed with reflection about the vertical axis acts on \(s=(q_{1x},q_{1y},q_{2x},q_{2y})\), and likewise on the momentum \(y\), by
\[
S=\operatorname{diag}(1,-1,-1,1),\qquad (s,y)\mapsto(Ss,Sy).
\]
For \(m_1=m_2\), \(F(Ss,Sy)=(SF_s,SF_y)\) and \(Ss_0=s_0\). This descends to an involution of the zero-angular-momentum rotation quotient.

The radial vector is \(s_0\). A rotation tangent vector is \((0,a,-b,0)\). The two mass-orthonormal shape tangent directions, both perpendicular to the radial/rotation directions, are
\[
e_+=(b,0,0,-a),\qquad Se_+=+e_+,
\]
\[
e_-=(0,b,a,0),\qquad Se_-=-e_-.
\]
The reduced second shape variation of \(U(s)=\sum m_i m_j/|a_{ij}\cdot s|\) in these directions is
\[
d_\pm=D^2U(s_0)[e_\pm,e_\pm].
\]
By the explicit pairwise Hessian identity
\[
D^2\left(\frac g{|r|}\right)[\Delta r,\Delta r]
=g\left(\frac{3(r\cdot\Delta r)^2}{|r|^5}
-\frac{|\Delta r|^2}{|r|^3}\right)
\]
and the three equal-side distances \(r_{12}=r_{13}=r_{23}=1/\rho_0\), one obtains
\[
\boxed{
d_- - d_+
=
\frac{3\rho_0(2M^2-M-1)}{M+2}
=
\frac{3\rho_0(2M+1)(M-1)}{M+2}>0.
}
\]

## Spectral assignment and forbidden terms

In either shape direction, the local \(2\times2\) linearized block from T16d is
\[
B(d)=\begin{pmatrix}-\nu_L&1\\d&\nu_L/2\end{pmatrix}.
\]
For fixed \(\nu_L<0\), the positive eigenvalue is strictly increasing in \(d\). Since \(d_->d_+\), the strong unstable eigenmode \(\xi_1\) and its paired stable eigenmode \(\xi_3\) are **odd** under \(S\), while the weak unstable mode \(\xi_2\) and paired stable mode \(\xi_4\) are **even**:
\[
S:(\xi_1,\xi_2,\xi_3,\xi_4)\mapsto(-\xi_1,+\xi_2,-\xi_3,+\xi_4).
\]
As the reduced nonlinear Newton field commutes with \(S\), each Taylor coefficient must respect these parities. But:
- \(\xi_2^2\) is even while \(\partial_{\xi_1}\) targets an odd vector component. Hence \(c_{\rm res}=0\).
- \(\xi_2\xi_3\) is odd while \(\partial_{\xi_4}\) targets an even component. Hence \(d_{\rm res}=0\).

These conclusions are independent of nonzero eigenvector rescalings and explain the T17e numerical observation \(|\kappa|<10^{-52}\) in the equal-first-two-masses family. The cancellation does **not** extend by symmetry alone to genuinely asymmetric mass triples (e.g. \((1,2,M_*)\)), where the reported nonzero values remain NUM pending interval validation.

## Audit

\`scripts/check_t21.py\` symbolically evaluates the three Newton potential Hessian contributions and proves their difference equals the boxed exact expression. The CI checks the derivation, **not** the independent correctness of the physical research context or literature originality.

This note corrects a preliminary hand-calculation typo \(3\rho_0(2M^2+5M-1)/(M+2)\): **that intermediate coefficient was false**; the boxed \(2M^2-M-1\) version is the independently symbolically checked one.
