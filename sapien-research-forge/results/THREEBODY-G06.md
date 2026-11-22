# G06 / T29 — Source-measure and endpoint corrections to three-body rare-event spectra

**Status:** Two exact mathematical theorems for explicitly specified local/toy models. Applying either as a global Newtonian three-body claim is **conditional** on source-measure and scattering regularity. Research lead Burak Çakay; no novelty priority asserted.

## T29a: source-measure exponent theorem

Let \(\Theta'_j=\lambda_j\Theta_j\), all \(\lambda_j>0\), with box exit \(|\Theta_j|=\delta_j>0\). Then
\[
\{T>t\}=\prod_j[-\delta_j e^{-\lambda_jt},\delta_j e^{-\lambda_jt}].
\]
Suppose \(f(\Theta)=C\prod_j|\Theta_j|^{\beta_j}(1+o(1))\) as \(\Theta\to0\), with uniform relative error, \(C>0,\ \beta_j>-1\). Integrating the product monomial over the exact box gives
\[
\boxed{P(T>t)\sim
C\prod_j\frac{2\delta_j^{\beta_j+1}}{\beta_j+1}
e^{-t\sum_j(\beta_j+1)\lambda_j}.}
\]
Thus, for \(r=e^{-v t}\) and \(v>0\), the source-dependent tail exponent is
\[
\boxed{\gamma_{\rm source}=\frac{\sum_j(\beta_j+1)\lambda_j}{v}.}
\]
Conditioned on \(T>t\), \(|\Theta_j|/(\delta_j e^{-\lambda_jt})\) has asymptotic CDF \(u^{\beta_j+1}\), independent across j. The remaining exit waits \(-\ln u_j/\lambda_j\) are therefore independent exponential variables with rates \((\beta_j+1)\lambda_j\), proving
\[
\boxed{P(\mathrm{face}=j\mid T>t)\to
\frac{(\beta_j+1)\lambda_j}{\sum_i(\beta_i+1)\lambda_i}.}
\]
T22/T23 are recovered by setting all \(\beta_j=0\), i.e. a smooth strictly positive source density at the unstable origin. A density vanishing on a mode can change the exponent *without changing the three-body masses or spectrum*. A singular lower-dimensional ensemble requires an additional tangent-rank argument and is **not covered** by this density formula.

For the numerical T23 spectrum \(m=(1,1,10)\), \(\lambda_j/v=(1.4226039399558574,0.75)\), the regular-density exponent is 2.1726039399558577, whereas \(f\propto\Theta_{\rm strong}^2\) gives
\(\gamma=3(1.4226039399558574)+0.75=5.017811819867572\) and strong-face probability 0.8505324577875872. These are local-model predictions, not astrophysical frequency claims.

## T29b: endpoint-dominance theorem in a conditional transfer toy model

Assume near an outgoing arc endpoint the weight obeys \(\mathfrak m(\psi)\asymp|\psi-\psi_b|^p\), \(p>0\), angular source density is regular, and the conditional close-approach probability is uniformly comparable to \(\min(1,(r/\mathfrak m(\psi))^\gamma)\). Then the endpoint integral is comparable to
\[
I(r)=\int_0^1\min(1,r^\gamma x^{-p\gamma})\,dx.
\]
Writing \(a=r^{1/p}\) yields the **exact identity**
\[
I(r)=a+r^\gamma\frac{1-a^{1-p\gamma}}{1-p\gamma},\quad p\gamma\ne1,
\]
and \(I(r)=r^\gamma[1+\log(1/r)/p]\) if \(p\gamma=1\). Hence
\[
\boxed{\gamma_{\rm effective}=\min(\gamma,1/p),\quad
p\gamma=1:\ I(r)\asymp r^\gamma\log(1/r).}
\]
This theorem does **not** assert that the Newtonian \(\mathfrak m(\psi)\) actually vanishes. The T27 certified pieces are interior subarcs; the endpoint sectors remain unverified. A positive uniform lower bound for \(\mathfrak m\) near the ends would *rule out* this anomaly.

## Independent Round20 archive audit

- ZIP supplied by Burak Çakay contains 231 entries, including code, an internally consistent pre-run prediction JSON and interval-output summaries. SHA-256 of archive: \`42d5b8c9746d3fc0636014f3f2858877b54520d3962700bfb687d961d1a7bf20\`.
- The prediction hash matches JSON serialization *after removing the hash field* using \`json.dumps(indent=1)\`. This confirms internal contents, **not** independent pre-run timestamp anchoring.
- Independently reran the original cocycle script: 19,136 pair/triple checks passed. Independently reran original M10 Newton trajectory driver for 250 cases: 250 events, 227 armed exits. The original 50k case run was not rerun.
- Independently enumerated all 40 degree-two shape-eigenvalue relations by exact algebra: precisely four relations, two each at \(k=(48\pm\sqrt{129})/75\).
- The supplied 141 marked-success interval pieces exactly cover the rational degree intervals \([1/256,93/32]\cup[124,23039/128]\) without gaps; 50 failed intermediate subdivision files also appear. All 141 named success files exist and their reported enclosures match the aggregate JSON.
- The weight integral interval arithmetic for the covered set independently checks out from its reported bounds:
  \(Z_w\in[23.2800,23.4938]\), \(Z_\ell\in[74.4647,77.5138]\), giving \(P_{\mathrm{short,covered}}\in[0.23096659,0.23983422]\).
- **Blocking upstream input:** the script \`r20_mweights.py\` and underlying \`r14_lohner.py\` require \`r10_manifold_certificate.json\`, **missing from the ZIP**. Thus the 141 interval calculations are still **CAP_U**, *not independently proof-checked*.
- The full original per-trajectory outputs \`r20_mc_M10.json\`, \`r20_mc_L.json\`, etc. were not included. Aggregate numerical outcome tables cannot be fully re-created from the ZIP alone.
- Pooled SE fields in \`r20_analyse.py\` are explicitly named \`se_if_independent\`. They combine overlapping windows of the *same run* and must not be interpreted as independent-trial confidence without joint/run-level bootstrap.
- The Arb-labelled point samples of the resonance coefficient \(c(\epsilon)\) for tiny \(\epsilon\) support \(c\propto\epsilon\), but do not by themselves rigorously bound \(c'(0)\) or certify nonzero throughout a full punctured mass interval. See P29 below.

## Near-term solution gates

1. **P28 (CAP replay):** provide \`r10_manifold_certificate.json\` plus its generator and verifier; independently rerun the 141 Lohner pieces with directed rounding, verify every terminal interval and zero gaps.
2. **P29 (resonant derivative):** prove \(0\notin\partial_\epsilon c_{\rm res}(0)\) by validated parameter-dependent eigenprojections and constrain Taylor remainder on \(|\epsilon|\le\epsilon_0\); thereby show/log-refute activation for *every* sufficiently small exchange asymmetry, not just sampled masses.
3. **P30 (physical source measure):** map a declared physical phase-space ensemble into the exact T17 eigenmode coordinates, certify its transverse rank and origin density/vanishing order; compare the appropriate T29a exponent on held-out runs.
4. **P31 (endpoint weight):** certify \(\inf\mathfrak m(\psi)>0\) and \(\inf F(\psi)>0\) on a defined full outgoing angular set, or prove a singular endpoint scaling and derive its impact via T29b. Without this, full-arc statistical exponents are conditional.

P15/P16 and a global non-tautological exact three-body transfer operator remain OPEN.

See \`registry/T20_R20_G06_AUDIT_P28_P31.json\` for machine-native definitions and acceptance/falsification scope.