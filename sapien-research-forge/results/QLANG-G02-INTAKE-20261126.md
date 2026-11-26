# Q-LANG / XENOGENESIS: finite semantics audit intake

Date: 2026-11-26. Status: SOURCE_REPORTED / INDEPENDENT_REPRODUCTION_PENDING. No novelty claim.

Source identities (summary only): `Q_LANG_REDTEAM_20261126.txt`, review draft, core G0.2, wire 0.1.1; `Q-LANG Manifest v0.1 — Formel Denetim Raporu.md`, 2026-11-25. Source documents were accessible to intake; original documents, private conversations and code are not published here. The latter report identifies canonical manifest SHA-256 `aec1304597bb040c1368c195f81783db2cca713d0b3c75320eecd69385e6a72f` (source-reported, not recomputed here).

## QL-005: finite symmetry check
For X={a,b} with transposition g, an equivariant singleton selector P cannot exist: gX=X implies gP(X)=P(X), but neither singleton is fixed. A nonempty set-valued selector P(X)=X is possible; an equivariant probability distribution on X must assign 1/2 to each atom. Source reports solver checks for these variants. Elementary mathematical observation, not a novelty claim; solver runs not independently repeated in this intake. Entropy of the distribution does not imply an external information lower bound.

## QL-G02: restricted semantics experiment
Reported profile: finite declared atoms A, actual domain D subset A (possibly empty), strict equality requiring both operands to denote, classical Boolean connectives, vacuous universal truth and false existential on empty domains. Reflexivity requires denotation; quantifier inference unsupported. Source reports 31/31 tests, 294 connective cases and 100 rejected mutations. These are reported measurements, not a fresh replay. Reproduce using original G0.2 archive and pinned environment before promotion.

## QL-006: evidence correction
Three SMT2 UNSAT results and a 512-relation finite scan are reported. Original AST, translator and independently checked proof terms remain missing or unverified. Backend entailment is not a checked proof. No proof status promotion. G0 remains OPEN; G1 locked.

## Reported limitations (not rerun)
- An expression depending on two variables can return UNKNOWN when 16 or 20 variables are declared; performance/coverage issue, not proof of unsoundness.
- A false evaluated assertion can coexist with a verified flag; distinguish evaluation success from assertion truth.
- An atom named `rat` reportedly collides with the AST parser.

No physical or metaphysical realization follows from these finite checks. No changes to P05 or X10. Human review required before merge.
