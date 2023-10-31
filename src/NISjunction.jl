U(ϵ, Δ) = 1/√(2)*√(1+√((ϵ^2-abs2(Δ))/ϵ^2))
V(ϵ, Δ) = 1/√(2)*√(1-√((ϵ^2-abs2(Δ))/ϵ^2))
aUV(U, V, Z) = U*V/(U^2+(U^2-V^2)*Z^2)
bUV(U, V, Z) = -(U^2-V^2)*(Z^2+im*Z)/(U^2+(U^2-V^2)*Z^2)
a(ϵ, Δ, Z) = aUV(U(ϵ, Δ), V(ϵ, Δ), Z)
b(ϵ, Δ, Z) = bUV(U(ϵ, Δ), V(ϵ, Δ), Z)

abs2a(ϵ, Δ, Z) = abs(ϵ) > abs(Δ) ? abs2(a(ϵ, Δ, Z)) : abs2(Δ)/(ϵ^2+(abs2(Δ)-ϵ^2)*(1+2Z^2)^2)
abs2b(ϵ, Δ, Z) = abs(ϵ) > abs(Δ) ? abs2(b(ϵ, Δ, Z)) : 1 - abs2a(ϵ, Δ, Z)

f0(ϵ, T) = T == 0 ? (ϵ < 0 ? 1 : 0) : 1/(1+exp(ϵ*ustrip(Unitful.q)/(ustrip(Unitful.k)*T)))
# If1(ϵ, Z, Δ) = 1-abs2(b(Z, ϵ, Δ))+abs2(a(Z, ϵ, Δ))
# If(ϵ, p; V=0, Z=0, Δ=1e-4*ustrip(Unitful.q), T=0) = (1-abs2(b(Z, ϵ, Δ))+abs2(a(Z, ϵ, Δ)))*(f0(ϵ-ustrip(Unitful.q)*V, T)-f0(ϵ, T))
function INIS(V; Z=1, A=1, N0=1e23, vF=100, S=1, T=0, Δ=1)
    e = ustrip(Unitful.q)
    kB = ustrip(Unitful.k)
    C = A*e*N0*vF*S
    If(ϵ, p) = C*(1-abs2b(ϵ, Δ, Z)+abs2a(ϵ, Δ, Z))*(f0(ϵ-V, T)-f0(ϵ, T))
    lb = min(0, V) - 3kB*T/e
    rb = max(0, V) + 3kB*T/e
    I = solve(IntegralProblem(If, lb, rb), QuadGKJL()).u
    return I
end

INIS(0.5; Z=50, T=0.01)

myΔ = 1e-5
Vs = (0:0.001:4)*myΔ
plot(Vs, (x->INIS(x; Z=10, Δ=myΔ)).(Vs))
plot(Vs, (x->INIS(x; Z=100, Δ=myΔ)).(Vs))





### zxy
