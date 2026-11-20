# G02: Finite transition representation for the reduced spatial three-body flow

**Status: OPEN_PROBLEM.** This is a research formulation, not a solved theorem.

**Research record:** Burak Çakay (@burakcakay), G02, in Sapien Research Forge. When this formulation is actually used, cite this precise file and the Git commit or accessed version.

## Problem

Let `m_i>0`, `ell=|J|>0`, `c=E ell^2` outside the relative-equilibrium critical values, and `X_c` a regular seven-dimensional translation/rotation-reduced, collision-free energy-angular momentum stratum with flow `Phi_t`.

For transverse event sections `Sigma_alpha`, define the exact transition only **implicitly** by

```math
P_{\alpha\beta}(z)=\Phi_{\tau_{\alpha\beta}(z)}(z),
\qquad g_\beta(P_{\alpha\beta}(z))=0.
```

Seek a finite, precisely delimited representation grammar `F` and event-internal interpolants `I_alpha` describing all required exact transitions **without including the Newton flow itself as a primitive**.

A proposed geometry-time transport, if the relevant rescaling conventions are held fixed, obeys

```math
(h_1,\phi_1,T_1)\star(h_2,\phi_2,T_2)
=(h_1+h_2,\phi_1+\phi_2,T_1+e^{3h_1/2}T_2).
```

The algebra of these reconstruction outputs does **not** itself provide a closed form of `P_alpha_beta`.

## Obstacles / falsification

1. Fix a non-tautological finite function grammar `F`.
2. Specify collision-free coverage, missed events, infinite-event accumulation and chart switching.
3. Include the known Lagrange equilateral Kepler relative periodic family as an exact special case.
4. Prove a construction in that grammar or give a rigorous impossibility result for that grammar.

No general closed-form flow-free operator is asserted. Do not cite this problem statement as a solution of the Newtonian three-body problem.
