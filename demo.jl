using ClassicalOrthogonalPolynomials, CairoMakie

T = ChebyshevT()
U = ChebyshevU()

U \ T

P = Legendre()
Q = Jacobi(1,0)

R = Q \ P
L = P \ Weighted(Q)

L*R

X = jacobimatrix(P)
I-X