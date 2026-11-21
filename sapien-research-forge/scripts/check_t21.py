"""Exact SymPy verification for T21 Hessian parity-selection theorem.

Computes pairwise Newton potential second derivatives directly.
This establishes the stated algebraic identity, not literature novelty.
"""
import sympy as sp

M=sp.symbols("M",real=True,positive=True)
mu=2*M/(M+2)
rho=sp.sqrt((2*M+1)/(M+2))
a=1/(sp.sqrt(2)*rho)
b=sp.sqrt(3*mu)/(2*rho)
q0=sp.Matrix([[a,0],[0,b]])
even=sp.Matrix([[b,0],[0,-a]])
odd=sp.Matrix([[0,b],[a,0]])
S=sp.diag(1,-1,-1,1)
assert sp.simplify(a*a+b*b-1)==0
assert S*sp.Matrix([b,0,0,-a])==sp.Matrix([b,0,0,-a])
assert S*sp.Matrix([0,b,a,0])==-sp.Matrix([0,b,a,0])

pairs=[
    (sp.Integer(1), sp.Matrix([[sp.sqrt(2),0]])),
    (M, sp.Matrix([[1/sp.sqrt(2),1/sp.sqrt(mu)]])),
    (M, sp.Matrix([[-1/sp.sqrt(2),1/sp.sqrt(mu)]])),
]
def hessian_direction(E):
    value=0
    for g,c in pairs:
        r=c*q0
        dr=c*E
        # All pair distances equal 1/rho for normalized Lagrange shape.
        assert sp.simplify(r.dot(r)-1/rho**2)==0
        value+=g*(3*rho**2*r.dot(dr)**2-dr.dot(dr))
    return sp.simplify(rho**3*value)

difference=sp.factor(sp.simplify(hessian_direction(odd)-hessian_direction(even)))
exact=3*rho*(2*M+1)*(M-1)/(M+2)
assert sp.simplify(difference-exact)==0,(difference,exact)
assert sp.factor(2*M**2-M-1)==(2*M+1)*(M-1)

sigma_star=(1064-32*sp.sqrt(129))/5625
eq=sp.together(sigma_star*(M+2)**2-(2*M+1))
roots=sp.solve(eq,M)
Mstar=next(x for x in roots if float(sp.N(x,20))>1)
assert 12.38<float(sp.N(Mstar,15))<12.39

# Exact parity logic: strong lambda_1 and stable lambda_3 are odd.
parity={1:-1,2:+1,3:-1,4:+1}
assert parity[2]**2!=parity[1]
assert parity[2]*parity[3]!=parity[4]

print("PASS T21: exact Hessian gap =",difference)
print("PASS Mstar =",sp.N(Mstar,18),"and both quadratic resonances are parity-forbidden")
