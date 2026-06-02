using MultivariateOrthogonalPolynomials, BlockArrays, CairoMakie, StaticArrays

ℓ = Int(block(B))-1
k = blockindex(B)
m = iseven(ℓ) ? k-isodd(k) : k-iseven(k)