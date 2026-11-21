# T12b / T17 — The invariants vanish on the analytic unstable parameterization

**Status: ANALYTIC_PROOF_CANDIDATE, conditional on the existence and analyticity of T10; independent human audit pending.**

Work with user-supplied equal-mass planar McGehee equations near collapsing Lagrange point \(L_-\). Let \(z(\theta)\) be analytic in an open neighborhood of \(\theta=0\in\mathbb R^2\), \(z(0)=L_-\), satisfying
\[
\Lambda(\theta\cdot\nabla_\theta)z(\theta)=F(z(\theta)),\qquad
\Lambda=\tfrac{\sqrt6}{4}(1+\sqrt{13})>0.
\]
Assume the exact differential identities
\[
F\cdot DG=-2\nu G,\qquad
F\cdot Dc=-(\nu/2)c,\qquad
F\cdot DH=\nu H
\]
and \(G(L_-)=c(L_-)=H(L_-)=0\), \(\nu(L_-)=-\sqrt6\).

**Proposition.** \(G\circ z=c\circ z=H\circ z\equiv0\) throughout the connected domain of analyticity where the parameterization identity is valid.

**Proof.** For \(K=G,c,H\), define the analytic scalar \(k(\theta)=K(z(\theta))\). The chain rule gives
\[
\Lambda(\theta\cdot\nabla_\theta)k(\theta)=\alpha_K(\theta)\,k(\theta),
\]
where \(\alpha_G(0)=2\sqrt6\), \(\alpha_c(0)=\sqrt6/2\), and \(\alpha_H(0)=-\sqrt6\).
If \(k\not\equiv0\), take the lowest nonzero homogeneous Taylor term \(k_n(\theta)\), degree \(n\in\mathbb N_{\ge1}\). Comparing homogeneous degree \(n\) yields
\[
n\Lambda k_n=\alpha_K(0)k_n.
\]
Therefore \(n=\alpha_K(0)/\Lambda\), but these three ratios are
\[
\frac{8}{1+\sqrt{13}}\in(1,2),\qquad
\frac{2}{1+\sqrt{13}}\in(0,1),\qquad
-\frac{4}{1+\sqrt{13}}<0.
\]
None is a positive integer. Contradiction. Thus each \(k\equiv0\) as a convergent Taylor series, hence on its connected analytic continuation domain. QED.

**Scope:** This eliminates an invariant-consistency gap in T12b, *assuming* T10's exact analytic parameterization exists with the stipulated invariance equation. It is not a validation of T10's purported interval coefficients, convergence radius or remainder.
