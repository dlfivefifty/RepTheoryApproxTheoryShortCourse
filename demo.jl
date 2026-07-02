# available on github
using ClassicalOrthogonalPolynomials, CairoMakie

T = ChebyshevT()
chebyshevt(5, 0.1) ≈ cos(5acos(0.1)) ≈ T[0.1,6]
plot(T[:,1:5])

U = ChebyshevU()

# raising operators
R = U \ T
# says T == U*R
# eg. (U[x,3]-U[x,1])/2 == T[x,3]
# eg. (U_2(x)-U_0(x))/2 == T[x,3]

U\diff(T)

P = Legendre()
Q = Jacobi(1,0)

R = Q \ P # raises just a in P^(a,b)
# We also have lowering operators between
# weighted Jacobi
W = Weighted(Q)
plot(W[:,4])
L = P \ W

# consider raising then lowering:
L*R

X = jacobimatrix(P)
I-X # == L*R

# L*R is the LU Factorisation of this polynomials weight modification
# with orthonormal it would sym tridiagonal
# and so LU becomes Cholesky factorisaytiomn
# Gautschi (1970) 
# showed computing cholesky of weight modification
# gives raising (change-of-basis/connection) matrices

# see also Gutleb, SO, Slevinsky (2024)
# for rational modifications