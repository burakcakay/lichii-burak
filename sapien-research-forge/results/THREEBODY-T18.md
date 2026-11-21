# T18 — A finite margin certificate for entire X10 balls

**Status: ANALYTIC_PROOF_CANDIDATE**, transparent algebraic enclosure; independent human review pending. This lemma proves a ball's containment in the X10 **sufficient escape cone**; it does not validate the upstream flow propagation producing the ball.

## Definitions

Fix a Newtonian *physical* planar Jacobi center \(\bar z=(\bar x,\bar X,\bar u,\bar W)\in\mathbb R^8\) and radius \(\epsilon>0\). Let
\[
B_\infty(\bar z,\epsilon)=\{z\in\mathbb R^8:\|z-\bar z\|_\infty\le\epsilon\}.
\]
Each of the four 2D blocks differs from its center by Euclidean norm at most \(\delta=\sqrt2\,\epsilon\). Fix \(b,V>0\) and put
\[
r_+=|\bar x|+\delta,\quad
R_-=|\bar X|-\delta,\quad
R_+=|\bar X|+\delta,\quad
d_-=R_- - b,\quad
u_+=|\bar u|+\delta.
\]
Define explicit conservative envelopes:
\[
e_+=\tfrac12u_+^2-\frac2{r_+},
\]
\[
L_-=|\bar x\wedge\bar u|-\delta(|\bar x|+|\bar u|)-\delta^2,
\]
\[
N_-=\bar X\cdot\bar W-\delta(|\bar X|+|\bar W|)-\delta^2,\qquad
s_-=\frac{N_-}{R_+}.
\]

**Finite certificate conditions** (only norms, dot/cross products, radicals, rational operations):
\[
\boxed{
|\bar x|-\delta>0,\quad b-r_+>0,\quad d_->0,\quad N_->0,
}
\]
\[
\boxed{
-\frac2b-e_+-\frac{2\sqrt b}{Vd_-^2}>0,\qquad
L_- -\frac{b^2}{Vd_-^2}>0,\qquad
s_-^2-V^2-\frac6{d_-}>0.
}
\]

## Theorem

If all seven quantities above are strictly positive, then
\[
\boxed{B_\infty(\bar z,\epsilon)\subset \Omega_{X10}(b,V)}
\]
in the physical, equal-mass planar Jacobi coordinates.

## Proof

Take any \((x,X,u,W)\in B_\infty(\bar z,\epsilon)\). Then
\[
|x|\le r_+<b,\quad |x|\ge|\bar x|-\delta>0,\quad
|X|\ge R_->b,\quad |X|\le R_+.
\]
Because \(e=\tfrac12|u|^2-2/|x|\) and \(|u|\le u_+\), \(e\le e_+\).
The bilinearity of the planar wedge yields
\[
|x\wedge u|\ge|\bar x\wedge\bar u|
-\delta(|\bar x|+|\bar u|)-\delta^2=L_-.
\]
The Cauchy inequality yields
\[
X\cdot W\ge
\bar X\cdot\bar W-\delta(|\bar X|+|\bar W|)-\delta^2=N_->0.
\]
Thus the radial velocity \(s=X\cdot W/|X|\ge N_-/R_+=s_->0\). Finally \(d=|X|-b\ge d_->0\). Monotonicity of the three remaining X10 margin expressions in their bounded arguments gives all their strict inequalities. QED.

**Use:** Apply T18 to each independently validated image ball at the exit of a T9–T15 interval step. This removes repeated pointwise checks inside the final X10 cone. It does **not** remove the need for validated integration, Levi-Civita chart transitions, and successful 141/141 coverage.

## Deterministic worked example

At \(\bar x=(0.3,0)\), \(\bar u=(0,1.5)\), \(\bar X=(8,0)\), \(\bar W=(3,0)\), \(b=0.41\), \(V=2\), and \(\epsilon=0.02\), the seven conservative conditions hold. The smallest of the seven margins is \(>0.0352\). This certifies all points of the full 8-dimensional cube, not merely its center.

**Limit:** In 3D, replace the planar scalar cross-product magnitude with the Euclidean norm of \(x\times u\). The same inequalities hold because the exterior product norm is submultiplicative.
