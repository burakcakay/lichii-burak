# X10 / T16 — A sufficient escape cone for the equal-mass Newtonian three-body problem

**Status: ANALYTIC_PROOF_CANDIDATE; explicit stand-alone derivation, independent human review pending.** This is a sufficient criterion, not a characterization of all escaping trajectories and not a newly established priority claim. Related Jacobi escape tests existed in the literature, including "Three-body problem—A test of escape valid even for very small mutual distances" (1984), https://www.sciencedirect.com/science/article/pii/0094576584900821 .

## Exact domain

Take equal positive masses m1=m2=m3=1, G=1, physical COM Jacobi coordinates
\[
x=r_2-r_1,\qquad X=r_3-\tfrac12(r_1+r_2),\qquad
u=\dot x,\qquad W=\dot X.
\]
The formulas hold in R^2 or R^3 (using the norm of vector angular momentum in 3D). Write \(g(y)=y/|y|^3,\ P=g(X+x/2),\ S=g(X-x/2)\). Newton's equations are
\[
\dot x=u,\quad\dot X=W,\quad
\dot u=-2x/|x|^3+S-P,\quad
\dot W=-\tfrac32(P+S).
\]
The scaled McGehee system in the chat has *different time and variable interpretations*; reconstruct physical (x,X,u,W) before applying these inequalities, or fix initial similarity scale rho=1.

Put
\[
r=|x|,\ R=|X|,\ e=\tfrac12|u|^2-2/r,\quad
L=x\wedge u,\ s=\dot R=X\cdot W/R.
\]
Choose \(b>r_0>0\), \(V>0\) with \(d=R_0-b>0\). Assume the exact strict conditions
\[
-\tfrac2b-e_0-\tfrac{2\sqrt b}{Vd^2}>0,\quad
|L_0|-\tfrac{b^2}{Vd^2}>0,\quad
s_0^2-V^2-\tfrac6d>0,\quad s_0>0.
\]

## Claim

For all physical \(t\ge0\), the solution exists and is collision-free, with
\[
\tfrac14 L_*^2<r(t)<b,\qquad s(t)>V,\qquad R(t)\ge R_0+Vt,\quad
L_*:=|L_0|-\tfrac{b^2}{Vd^2}>0.
\]
The outer velocity converges to \(W_\infty\), \(|W_\infty|\ge V\), and the outer asymptotic Kepler energy
\[
e_\infty=\tfrac12|W_\infty|^2\ge\tfrac12V^2>0.
\]
In the COM frame, \(r_3=(2/3)X\), so \(|r_3|\to\infty\). The pair \(1,2\) stays bounded without collision.

## Proof

**(1) Differential tidal bounds.** For any \(y\ne0\),
\[
Dg(y)=|y|^{-3}I-3|y|^{-5}yy^\mathsf T,
\quad\|Dg(y)\|_2=2|y|^{-3}.
\]
Thus, as long as \(r<b<R\),
\[
|S-P|\le 2r/(R-r/2)^3\le 2r/(R-b)^3,
\qquad |\dot W|\le3/(R-b)^2.
\]

**(2) Bootstrap radial motion.** On any interval with \(s>V\) and \(r<b\),
\[
\ddot R=(|W|^2-s^2)/R+(X/R)\cdot\dot W
\ge -3/(R-b)^2.
\]
Multiplying by \(s=\dot R>0\) and integrating gives
\[
s(t)^2\ge s_0^2-6/d+6/(R(t)-b)>V^2.
\]
It follows that \(R(t)\ge R_0+Vt\).

**(3) Bootstrap negative binary energy.** On the interval where \(e<0\),
\[
|u|^2<4/r,\quad
|\dot e|=|u\cdot(S-P)|\le4\sqrt b/(R-b)^3
\le4\sqrt b/(d+Vt)^3.
\]
Therefore \(|e(t)-e_0|\le2\sqrt b/(Vd^2)\). The stated margin implies \(e(t)<-2/b\); because \(e\ge-2/r\), this implies \(r(t)<b\). The bootstrap condition reproduces itself.

**(4) Bootstrap binary angular momentum.** Regardless of planar/spatial orientation,
\[
|\dot L|\le |x|\,|S-P|\le2b^2/(d+Vt)^3.
\]
Consequently \(|L(t)-L_0|\le b^2/(Vd^2)\), so \(|L(t)|\ge L_*>0\).
Since \(e<0\) also gives \(|u|^2<4/r\), Cauchy yields
\[
L_*^2\le |L|^2\le r^2|u|^2<4r
\quad\Rightarrow\quad r>L_*^2/4.
\]
The other two pairwise distances exceed \(R-r/2>R-b/2>0\).

**(5) Global existence and asymptotic escape.** All forces remain bounded on every finite interval, thanks to \(r>L_*^2/4\), \(R-r/2>0\), and \(\dot W\) being bounded. Thus the locally smooth Newton solution continues for all positive times. As
\[
\int_0^\infty|\dot W|\,dt
\le\int_0^\infty 3(d+Vt)^{-2}dt=3/(Vd),
\]
\(W(t)\) converges to \(W_\infty\). Since \(s(t)>V\), \(|W_\infty|\ge V\). Moreover \(R\to\infty\); the outer monopole potential tends to zero, yielding \(e_\infty=|W_\infty|^2/2\ge V^2/2>0\). QED.

**Strictness:** The proof works with exact strict inequalities and requires no finite-time numerical integration. It proves the *escape implication*, not that any T9–T15 interval box actually lies inside the cone. Those claims are still CAP_U until separately reproduced.

## Example / independent numerical crosscheck

A physical starting state \(x=(0.30,0)\), \(u=(0,1.5)\), \(X=(8,0)\), \(W=(3,0)\), \(b=0.41\), \(V=2\) satisfies the five margins approximately \((0.11,\ 7.59,\ 0.65250,\ 0.44854,\ 4.20949)\). A DOP853 numerical check with 1e-11 requested tolerances to \(t=100\) showed no collision, \(R(100)\approx296.52\), \(\min r\approx0.06090\), and \(\min s\approx2.87578\). **This numerical check is corroboration, not a computer-assisted proof.**
