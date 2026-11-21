# G05 / T20 — Spectral exit flux, conditional close-approach tails, and two corrections

**Status:** the local linear exit-flux formula is an **exact elementary theorem**. Its use to predict global Newtonian three-body statistics is **conditional**. This record **does not assert novelty**; related McGehee/collision and three-body statistical-flux literature exists.

**Provenance:** Research program directed by Burak Çakay (@burakcakay), based on T17–T19 claims supplied in conversation 2026-11-22, and the linked G03/G04 mathematics. Original \`r17_*\`, \`r18_*\`, \`r19_*\` calculation source was *not independently rerun*.

## 1. Exact local exit-flux theorem

Assume \(\Lambda_+,\Lambda_-,\delta>0\), \(\Theta'(t)=\operatorname{diag}(\Lambda_+,\Lambda_-)\Theta(t)\), and exit section \(|\Theta|=\delta\). Parameterize each initial point by its forward exit time \(T\) and exit angle \(\psi\):
\[
\Theta_0(T,\psi)=\delta(e^{-\Lambda_+T}\cos\psi,e^{-\Lambda_-T}\sin\psi).
\]
The change-of-variables determinant is
\[
\boxed{\det\frac{\partial(\Theta_{01},\Theta_{02})}{\partial(T,\psi)}
=-\delta^2e^{-(\Lambda_++\Lambda_-)T}
[\Lambda_+\cos^2\psi+\Lambda_-\sin^2\psi].}
\]
If the initial density \(f\) is continuous and strictly positive at the origin, then conditional on \(T>T_0\), as \(T_0\to\infty\), the exit-angle density is
\[
\boxed{p_{\rm exit}(\psi)=
\frac{\Lambda_+\cos^2\psi+\Lambda_-\sin^2\psi}
{\pi(\Lambda_++\Lambda_-)}.}
\]
This follows by inserting \(f(\Theta_0)\to f(0)\) uniformly in \(\psi\) and using
\(\int_0^{2\pi}(\Lambda_+\cos^2\psi+\Lambda_-\sin^2\psi)\,d\psi=\pi(\Lambda_++\Lambda_-)\).

For a *square* exit \(\max(|\Theta_1|,|\Theta_2|)=\delta\), the asymptotic probabilities of exiting through either pair of opposite faces are exactly
\[
\boxed{P_+=\Lambda_+/(\Lambda_++\Lambda_-),\qquad
P_-=\Lambda_-/(\Lambda_++\Lambda_-).}
\]
The face flux density is \(\delta\Lambda_i e^{-(\Lambda_++\Lambda_-)T}\,dT\,d(\text{other exit coordinate})\), hence independent of the other face coordinate to leading order.

In the T17 unequal-mass Lagrange spectrum with \(m=(1,1,10)\),
\[
\gamma=\frac{\Lambda_++\Lambda_-}{|\nu_L|}=2.1726039399558577,\quad
(P_+,P_-)=(0.6547921200882851,0.3452078799117148).
\]

## 2. Conditional global close-approach and ejection laws

Fix a collision-free outgoing branch on which a uniform return/exit map exists and, for some \(\kappa>0\), \(v=|\nu_L|\),
\[
\rho_{\min}=\kappa e^{-vT}\mathfrak m(\psi)(1+o(1)),
\qquad \mathfrak m(\psi)>0,\qquad \gamma=(\Lambda_++\Lambda_-)/v.
\]
Assume all errors are uniform on a compact angular window, and the initial transverse density is smooth and positive. The local flux identity implies the **conditional** asymptotics
\[
\boxed{P(\rho_{\min}<r,\psi\in W)\sim
C\,r^\gamma\int_W
a(\psi)\mathfrak m(\psi)^{-\gamma}\,d\psi,}
\quad a(\psi)=\Lambda_+\cos^2\psi+\Lambda_-\sin^2\psi.
\]
For a certified hyperbolic-escape branch on which \(F(\psi)=\rho_{\min}e_\infty>0\) is bounded and uniform in the limit,
\[
\boxed{P(e_\infty>E,\psi\in W)\sim
C\,E^{-\gamma}\int_W
a(\psi)[F(\psi)/\mathfrak m(\psi)]^\gamma\,d\psi.}
\]
If \(e_\infty=v_\infty^2/2\), the corresponding *speed tail* has exponent \(2\gamma\).

These are *not proven for the full three-body problem*, since uniform regularity of \(\mathfrak m,F\), branch continuation, energy fibers and endpoint exclusions remain open. For \(\Lambda_+=\Lambda_-\), the factor \(a(\psi)\) is constant and the T18a/b equal-mass weights are recovered.

## 3. Exact T17b chart-transition correction

The Hermitian planar shape matrix is \(H_{ij}=\bar z_i z_j\). On \(U_k=\{H_{kk}>0\}\), gauge \(z_k^{(k)}>0\) is \(z_j^{(k)}=H_{kj}/\sqrt{H_{kk}}\). Thus on overlaps,
\[
\boxed{z^{(\ell)}=\frac{\overline{H_{k\ell}}}{|H_{k\ell}|}\,z^{(k)},}
\]
**not** \(H_{k\ell}z^{(k)}/|H_{k\ell}|\). Witness \(z=(1,i)\) has \(H_{12}=i\), \(z^{(1)}=(1,i)\), \(z^{(2)}=(-i,1)\). This matters for exact symplectic/orientation-compatible chart gluing; the T17b orbit-space bijection remains valid.

## 4. Exact T18a finite-radius density counterexample

T18a.d stated \(P(\psi\in d\psi\mid|\Theta_0|<R)=d\psi/(2\pi)\) from only smooth positivity of initial density. This equality is **false at fixed \(R\)**. On \(|\Theta|<R\) take \(f(\Theta)=1+a\Theta_1\), with \(|a|R<1\). Direct polar integration gives
\[
\boxed{p_R(\psi)=\frac1{2\pi}\left(1+\frac{2aR}{3}\cos\psi\right).}
\]
The corrected universal statement is \(p_R(\psi)\to1/(2\pi)\) as \(R\to0\). Exact finite-\(R\) uniformity requires a special radial-symmetry condition or an equivalent angular cancellation.

## 5. Angular-momentum double scaling (not yet a tail theorem)

The exact McGehee angular momentum is \(c=J/\sqrt\rho\). Thus staying in a bounded transverse angular-momentum chart near \(\rho\to0\) requires \(|J|=O(\sqrt\rho)\). For a closest-approach cutoff \(r\), the natural crossover variable is
\[
\boxed{\zeta_L=J/\sqrt r.}
\]
If a regular planar energy-section ensemble has a genuinely independent transverse \(J\) coordinate, with uniform outgoing controls, its additional small window width is \(O(r^{1/2})\), suggesting \(\gamma_{\rm planar,freeL}=\gamma+1/2\). A spatial 3-component vector would instead suggest \(\gamma+3/2\) **only after proving rank 3 and uniformity**.

The T19 free-L numerical pilot \(1.05\pm0.09\) versus proposed \(2.8028\) is a substantive discrepancy, *not* proof of the latter. Likewise the \(m=(1,1,10)\) pilot \(1.38\pm0.11\) versus \(2.1726\) remains unresolved. Larger sample size alone cannot fix an unverified source measure or incorrect asymptotic limit.

## 6. Status inheritance from supplied T17–T19

- T17a: high-precision **NUM**, Euler time factor \(+\nu_E/4\), not Lagrange \(\nu_L/4\).
- T17b: **SYM** for Hermitian rank-one orbit-space map, with gauge transition corrected above.
- T17c: **SYM+NUM** conjugacy of G03 quartic flow and planar McGehee dynamics.
- T17d: **formal resonant classification**; analytic convergence and quantitative exit map remain open.
- T17e: symplectic coefficient identity **PP_U**, numerical nonzero \(\kappa\) only, interval proof open; equal-pair mass symmetry cancellation is **NUM**.
- T18a/b: proof candidates subject to measure/branch and uniformity assumptions.
- T18c: free angular momentum/spatial exponent remains **COND**.
- T18d/e: calculated weights and MC evidence are **NUM**.
- T19: \(h\neq0\) exponent \(2.22\pm0.18\) (weakly consistent), free \(L\) exponent \(1.05\pm0.09\), \(m=(1,1,10)\) exponent \(1.38\pm0.11\) (both do not confirm predicted scaling).
- P14 coverage reported as **CAP_U** 141/141 at \(\theta_m=.003\), plus 29/29 newly reported short window; \(353.3672^\circ\) of 360 reported, still unaudited. P15/P16 remain **OPEN**.
- Larger T19 runs were *not supplied or evaluated* for this record.

## 7. Related work to compare before asserting novelty

- [Three-body scattering near triple collision or expansion (MNRAS, 1991)](https://academic.oup.com/mnras/article/250/3/555/998448).
- [Testing the flux-based statistical prediction of the three-body problem (MNRAS, 2021)](https://academic.oup.com/mnras/article/506/1/692/6297275).

The local Jacobian derivation above is elementary. A genuinely publishable contribution would require a **global verified branch-wise law**, auditable uniform bounds and independent observational/Monte-Carlo falsification tests.
