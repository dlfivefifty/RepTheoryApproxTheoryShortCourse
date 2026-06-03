using MultivariateOrthogonalPolynomials, BlockArrays, CairoMakie, StaticArrays, LaTeXStrings, ClassicalOrthogonalPolynomials
using BlockArrays: block, blockindex

yₘ = m -> splat((x,y) -> m ≥ 0 ? (x+im*y)^abs(m) : (x-im*y)^abs(m))


𝐲ₘ = m -> function(𝐱)
    r = norm(𝐱)
    θ = atan(reverse(𝐱)...)
    r^(abs(m)-1)*exp(im*m*θ)*(𝐞_r(𝐱) + (1-2signbit(m))*im*𝐞_θ(𝐱))
end

Σ = splat((x,y) -> [x^2-y^2 2x*y; 2x*y y^2-x^2])


𝐞_r = splat((x,y) -> SVector(x,y) / sqrt(x^2+y^2))
𝐞_θ = splat((x,y) -> SVector(-y,x) / sqrt(x^2+y^2))

zₘ = function (m,λ,j,(x,y))
    r² = x^2+y^2
    jacobip(j,λ,abs(m),2r²-1) * yₘ(m)(SVector(x,y))
end


𝐳ₘ¹ = function (m,λ,j,(x,y))
    r² = x^2+y^2
    jacobip(j,λ,abs(m)-1,2r²-1) * 𝐲ₘ(m)(SVector(x,y))
end
𝐳ₘ² = function (m,λ,j,(x,y))
    r² = x^2+y^2
    jacobip(j-1,λ,abs(m)+1,2r²-1) * Σ((x,y)) * 𝐲ₘ(m)(SVector(x,y))
end

P = (n,a,b,x) -> jacobip(n,a,b,x)
𝐩 = (n,a,b,α,t) -> if α == 1
    P(n, a, b, 2t-1) * [1,1]
else
    t*P(n-1, a, b+2, 2t-1) * [1,-1]
end

𝐏 = (n,a,b,t) -> n == 0 ? 𝐩(n,a,b,1,t) : [𝐩(n,a,b,2,t) 𝐩(n,a,b,1,t)]

𝔯 = (n,b) -> 10n^2 + n*(13b+19) + 4*(b+1)*(b+2)
𝐯 = (n,b,t) -> [-2t*(n+1)*P(n,0,b+2,2t-1) + (b+1)*(1-t)*P(n,1,b+1,2t-1); (b+1)*P(n,1,b+1,2t-1)]
𝐪 = (n,b,α,t) -> if α == 1
        n < 0 && return zero(t)
        (2n+b+1)/((n+1)*𝔯(n,b)) * ((2n+b+2)*(2n+b+3)*[(1-t)*P(n,1,b,2t-1); P(n,1,b,2t-1)] - (n+b+2) * 𝐯(n,b,t))
    else
        n < 1 && return zero(t)
        -(𝐯(n-1,b,t)/n + 𝐪(n,b,1,t))
    end
pᴺ = function(m,j,α,𝐱)
    if iszero(m)
        (x,y) = 𝐱
        r² = x^2+y^2
        α == 1 ? P(j-1,0,1,2r²-1)*[x,y] : P(j-1,1,1,2r²-1)*[-y,x]
    else
        𝒫(m)(t -> 𝐪(j,abs(m)-1,α,t),𝐱)
    end
end
𝐧 = (m,j,α,(x,y)) -> [1-y^2 x*y; x*y 1-x^2]*pᴺ(m,j,α,SVector(x,y))
𝐧⁺ = (m,j,𝐱) -> m == 0 ? 2𝐧(m,j,1,𝐱) : 𝐧(m,j,1,𝐱) + 𝐧(m,j,2,𝐱)
𝐧⁻ = function(m,j,𝐱)
    if m == 0
        2𝐧(0,j,2,𝐱)
    elseif j == 0
        2𝐧(m,0,1,𝐱)
    else
        𝐧(m,j,1,𝐱) - 𝐧(m,j,2,𝐱)
    end
end

function blockindex2mj(B::BlockIndex)
    ℓ = Int(block(B))-1
    k = blockindex(B)
    m = iseven(ℓ) ? k-isodd(k) : k-iseven(k)
    (isodd(k+ℓ) ? 1 : -1) * m, (ℓ-m) ÷ 2
end
Z = Zernike()
for B in Block(5)[1:5]
    m,j = blockindex2mj(B)
    fig = Figure(size = (800, 800)); ax = Axis(fig[1,1], aspect = 1.0)
    contourf!(ax, Z[:,B]); text!(0.6, 0.8;  text=LaTeXString((m == 0 ? "" : (m > 0 ? "Re " : "Im " )) *  "z\${}_{$(abs(m))$j}\$"), fontsize=55); fig
    save("figures/z_$m$j.png", fig)
end


W = Weighted(Zernike(1))
for B in Block(5)[1:5]
    m,j = blockindex2mj(B)
    fig = Figure(size = (800, 800)); ax = Axis(fig[1,1], aspect = 1.0)
    contourf!(ax, W[:,B]); text!(0.55, 0.82;  text=LaTeXString((m == 0 ? "" : (m > 0 ? "Re " : "Im " )) *  "w\${}_{$(abs(m))$j}\$"), fontsize=55); fig
    save("figures/w_$m$j.png", fig)
end


m = 2; j = 1;
# contourf(Z*real(Z\Weighted(DiskBubble(m))[:,j+1]))
fig = Figure(;size=(800,800))
ax = Axis(fig[1,1]; aspect=1) 
streamplot!(ax, 𝐱 -> norm(𝐱) ≤ 1 ? real(Weighted(DiskNormal(m))[𝐱,4]+Weighted(DiskNormal(m))[𝐱,5]) : NaN*𝐱, -1..1, -1..1, arrow_size=20, gridsize=(35,35), linewidth=4)
text!(ax, 0.55, 0.82;  text=LaTeXString("Re∇w\${}_{$m$j}\$"), fontsize=50)
fig
save("figures/∇w_$m$j.png", fig)

fig = Figure(;size=(800,800))
ax = Axis(fig[1,1]; aspect=1) 
streamplot!(ax, 𝐱 -> norm(𝐱) ≤ 1 ? imag(Weighted(DiskNormal(m))[𝐱,4]+Weighted(DiskNormal(m))[𝐱,5]) : NaN*𝐱, -1..1, -1..1, arrow_size=20, gridsize=(35,35), linewidth=4)
text!(ax, 0.55, 0.82;  text=LaTeXString("Im∇w\${}_{$m$j}\$"), fontsize=50)
fig
save("figures/∇w_-$m$j.png", fig)



fig = Figure(;size=(800,800))
ax = Axis(fig[1,1]; aspect=1) 
streamplot!(ax, 𝐱 -> norm(𝐱) ≤ 1 ? SVector{2}(𝐧(0,4,2,𝐱)+𝐧(0,3,1,𝐱)) : NaN*𝐱, -1..1, -1..1, arrow_size=20, gridsize=(35,35), linewidth=4)
fig
save("figures/equivar.png", fig)


fig = Figure(;size=(800,800))
ax = Axis(fig[1,1]; aspect=1) 
m,j = 3,2;
streamplot!(ax, 𝐱 -> norm(𝐱) ≤ 1 ? SVector{2}(real(𝐳ₘ¹(m,0,j,𝐱))) : NaN*𝐱, -1..1, -1..1, arrow_size=20, gridsize=(35,35), linewidth=4)
text!(ax, 0.55, 0.82;  text=LaTeXString("Re 𝐳\${}^1_{$m$j}\$"), fontsize=50)
fig
save("figures/re𝐳_$m$j^1.png", fig)


fig = Figure(;size=(800,800))
ax = Axis(fig[1,1]; aspect=1) 
m,j = 3,2;
streamplot!(ax, 𝐱 -> norm(𝐱) ≤ 1 ? SVector{2}(imag(𝐳ₘ¹(m,0,j,𝐱))) : NaN*𝐱, -1..1, -1..1, arrow_size=20, gridsize=(35,35), linewidth=4)
text!(ax, 0.55, 0.82;  text=LaTeXString("Im 𝐳\${}^1_{$m$j}\$"), fontsize=50)
fig
save("figures/im𝐳_$m$j^1.png", fig)


fig = Figure(;size=(800,800))
ax = Axis(fig[1,1]; aspect=1) 
m,j = 3,2;
streamplot!(ax, 𝐱 -> norm(𝐱) ≤ 1 ? SVector{2}(real(𝐳ₘ²(m,0,j,𝐱))) : NaN*𝐱, -1..1, -1..1, arrow_size=20, gridsize=(35,35), linewidth=4)
text!(ax, 0.55, 0.82;  text=LaTeXString("Re 𝐳\${}^2_{$m$j}\$"), fontsize=50)
fig
save("figures/re𝐳_$m$j^2.png", fig)


fig = Figure(;size=(800,800))
ax = Axis(fig[1,1]; aspect=1) 
m,j = 3,2;
streamplot!(ax, 𝐱 -> norm(𝐱) ≤ 1 ? SVector{2}(imag(𝐳ₘ²(m,0,j,𝐱))) : NaN*𝐱, -1..1, -1..1, arrow_size=20, gridsize=(35,35), linewidth=4)
text!(ax, 0.55, 0.82;  text=LaTeXString("Im 𝐳\${}^2_{$m$j}\$"), fontsize=50)
fig
save("figures/im𝐳_$m$j^2.png", fig)


t = range(0,1,100)
n = 0
fig = lines(t, first.(𝐪.(n,0,1, t)), linewidth=4)
lines!(t, last.(𝐪.(n,0,1, t)), linewidth=4)
text!(0.2, 0.82;  text=LaTeXString("\$𝐪_{$n}^1\$"), fontsize=50)
ylims!(-2,2)
fig

n = 3
fig = lines(t, first.(𝐪.(n,0,1, t)), linewidth=4)
lines!(t, last.(𝐪.(n,0,1, t)), linewidth=4)
text!(0.2, 0.82;  text=LaTeXString("\$𝐪_{$n}^1\$"), fontsize=50)
ylims!(-2,2)
fig

for n = 0:3, ν = 1:min(n+1,2)
    fig = lines(t, first.(𝐪.(n,0,ν, t)), linewidth=4)
    lines!(t, last.(𝐪.(n,0,ν, t)), linewidth=4)
    text!(0.2, 0.82;  text=LaTeXString("\$𝐪_{$n}^{(0),$ν}\$"), fontsize=50)
    ylims!(-2,2)
    save("figures/𝐪_$n^$ν.png", fig)
end