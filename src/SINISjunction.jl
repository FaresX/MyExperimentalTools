using Integrals
using LinearSolve
using LinearAlgebra
using IterativeSolvers
using Unitful
using GLMakie

const e = ustrip(Unitful.q)
const kB = ustrip(Unitful.k)

Ain(E; Δ=0.0001, Z=0) = Δ^2 / (E^2 + (1 + 2Z^2)^2 * (Δ^2 - E^2))
Bin(E; Δ=0.0001, Z=0) = 4Z^2 * (1 + Z^2) * (Δ^2 - E^2) / (E^2 + (1 + 2Z^2)^2 * (Δ^2 - E^2))
Tin(E) = 0

Aout(E; Δ=0.0001, Z=0) = Δ^2 / (E + (1 + 2Z^2) * √(E^2 - Δ^2))^2
Bout(E; Δ=0.0001, Z=0) = 4Z^2 * (1 + Z^2) * (E^2 - Δ^2) / (E + (1 + 2Z^2) * √(E^2 - Δ^2))^2
Tout(E; Δ=0.0001, Z=0) = 2(E^2 - Δ^2 + E * (1 + 2Z^2) * √(E^2 - Δ^2)) / (E + (1 + 2Z^2) * √(E^2 - Δ^2))^2

A(E; Δ=0.0001, Z=0) = abs(E) < Δ ? Ain(E; Δ=Δ, Z=Z) : Aout(E; Δ=Δ, Z=Z)
A(n, V, E; Δ=0.0001, Z=0) = A(E + n * e * V; Δ=Δ, Z=Z)
B(E; Δ=0.0001, Z=0) = abs(E) < Δ ? Bin(E; Δ=Δ, Z=Z) : Bout(E; Δ=Δ, Z=Z)
B(n, V, E; Δ=0.0001, Z=0) = B(E + n * e * V; Δ=Δ, Z=Z)
T(E; Δ=0.0001, Z=0) = abs(E) < Δ ? Tin(E) : Tout(E; Δ=Δ, Z=Z)
T(n, V, E; Δ=0.0001, Z=0) = T(E + n * e * V; Δ=Δ, Z=Z)
f0(E; μ=0, t=0) = 1 / (1 + exp((E - μ) * e / (kB * t)))
f(n, V, E; μ=0, t=0) = E - μ == t == 0 ? 1.0 : f0(E + n * e * V; μ=μ, t=t)
δ(i, j) = float(i == j)

function Q(V, E, i, j; Δ=0.0001, Z=0, n=10)
    return if i <= 2n + 1
        A(n + 1 - i, V, E; Δ=Δ, Z=Z) * δ(j, i + 1) - B(n + 1 - i, V, E; Δ=Δ, Z=Z) * δ(j, 4n + 5 - i)
    elseif i == 2n + 2
        δ(i, j)
    elseif 2n + 3 <= i <= 4n + 3
        A(3n + 3 - i, V, -E; Δ=Δ, Z=Z) * δ(j, i + 1) - B(3n + 3 - i, V, -E; Δ=Δ, Z=Z) * δ(j, 4n + 5 - i)
    elseif i == 4n + 4
        δ(i, j)
    end
end

function P(V, E, i; Δ=0.0001, Z=0, μ=0, t=0, n=10)
    return if i < 2n + 2
        B(n + 1 - i, V, E; Δ=Δ, Z=Z) + T(n + 1 - i, V, E; Δ=Δ, Z=Z) * f(n + 1 - i, V, E; μ=μ, t=t)
    elseif i == 2n + 2
        0.0
    elseif 2n + 3 <= i <= 4n + 3
        B(3n + 3 - i, V, -E; Δ=Δ, Z=Z) + T(3n + 3 - i, V, -E; Δ=Δ, Z=Z) * f(3n + 3 - i, V, -E; μ=μ, t=t)
    elseif i == 4n + 4
        0.0
    end
end

function G0(V, E, i; μ=0, t=0, n=10)
    return if i <= 2n + 2
        f(n + 1 - i, V, E; μ=μ, t=t)
    # elseif i == 2n + 2
    #     1.0
    # elseif 2n + 3 <= i <= 4n + 3
    else
        f(3n + 3 - i, V, -E; μ=μ, t=t)
    # elseif i == 4n + 4
    #     1.0
    end
end

cut(x) = 0 <= x <= 1 ? x : x < 0 ? 0.0 : 1.0

V = 1e-3
Δ = 1e-4
E = 1e-3
Z = 1
μ = 0
t = 0
n = 1
η = 0.01
A(1e-4)
P(1e-3, 1e-4, 1)
T(1e-3)
G0m = [G0(V, E, i; μ=0, t=0, n=n) for i in 1:4n+4]
Qm = [Q(V, E, i, j; Δ=Δ, Z=Z, n=n) for i in 1:4n+4, j in 1:4n+4]
Pm = [P(V, E, i; Δ=Δ, Z=Z, μ=μ, t=t, n=n) for i in 1:4n+4]
ΔG = cut.(Qm * G0m + Pm) - G0m
Gn = G0m + ΔG * η
actual_tol = 1
while true
    Go = copy(Gn)
    ΔG = cut.(Qm * Gn + Pm) - Gn
    Gn += ΔG * η
    actual_tol = sum(abs2.(Gn .- Go))
    actual_tol < 1e-30 && break
end
Gn
actual_tol

function gs(V, E; Δ=0.0001, Z=0, μ=0, t=0, n=10, iters=1e9, tol=1e-40)
    G0m = [G0(V, E, i; μ=μ, t=t, n=n) for i in 1:4n+4]
    Qm = [Q(V, E, i, j; Δ=Δ, Z=Z, n=n) for i in 1:4n+4, j in 1:4n+4]
    Pm = [P(V, E, i; Δ=Δ, Z=Z, μ=μ, t=t, n=n) for i in 1:4n+4]
    ΔG = cut.(Qm * G0m + Pm) - G0m
    Gn = G0m + ΔG * η
    actual_tol = 1
    for _ in 1:iters
        Go = copy(Gn)
        ΔG = cut.(Qm * Gn + Pm) - Gn
        Gn += ΔG * η
        actual_tol = sum(abs2.(Gn .- Go))
        actual_tol < tol && break
    end
    return Gn
end
go(E; Δ=0.0001, Z=0, μ=0, t=0, n=10, iters=1e9, tol=1e-40) = gd(V, E; Δ=Δ, Z=Z, μ=μ, t=t, n=n, iters=iters, tol=tol)[n+1]
function gd(V, E; Δ=0.0001, Z=0, μ=0, t=0, n=10, iters=1e9, tol=1e-40)
    Gs = gs(V, E; Δ=Δ, Z=Z, μ=μ, t=t, n=n, iters=iters, tol=tol)
    Gs[n+1] + Gs[3n+4]
end

gd(1e-3, 1.001e-4)

function INIS(V; Z=0, A=1, N0=1e23, vF=100, S=1, t=0, Δ=0.0001, μ=0, n=10, iters=1e9, tol=1e-40)
    e = ustrip(Unitful.q)
    kB = ustrip(Unitful.k)
    C = 2A * e * N0 * vF * S
    If(E, p) = C * (gd(V, E; Δ=Δ, Z=Z, μ=μ, t=t, n=n, iters=iters, tol=tol) - 1)
    lb = min(0, V) - 3kB * t / e
    rb = max(0, V) + 3kB * t / e
    I = solve(IntegralProblem(If, lb, rb), QuadGKJL()).u
    return I
end

INIS(1e-5; n=1)
V = -2e-4:1e-5:2e-4
Is = [INIS(v; n=1) for v in V]
lines(Is, V)



gsp(E; Z=1) = ((1 - A(E; Z=Z)) * B(E; Z=Z) - B(E; Z=Z)^2) / ((1 - A(E; Z=Z))^2 - B(E; Z=Z)^2)
gsp(1e-5)
(1 - A(1e-5))^2 - B(1e-5)^2
B(1e-3)

ft(E; Z=1) = B(E; Z=1) * B(-E; Z=1) / ((1 - A(E; Z=1)) * (1 - A(-E; Z=1)))
ft(1.05e-4)
A(1e-6)

1 - A(1E-3) ≈ B(1e-3)



#########################################
function g(E; Δ=0.0001, Z=0, μ=0, t=0)
    return if abs(E) <= Δ
        1 / 2
    else
        (1 - Aout(-E; Δ=Δ, Z=Z)) / ((1 - Aout(E; Δ=Δ, Z=Z)) * (1 - Aout(-E; Δ=Δ, Z=Z)) - Bout(E; Δ=Δ, Z=Z) * Bout(-E; Δ=Δ, Z=Z)) *
        (Bout(E; Δ=Δ, Z=Z) + Tout(E; Δ=Δ, Z=Z) * f0(E; μ=μ, t=t) - (Bout(E; Δ=Δ, Z=Z) * (Bout(-E; Δ=Δ, Z=Z) + Tout(-E; Δ=Δ, Z=Z) * f0(-E; μ=μ, t=t))) / (1 - Aout(-E; Δ=Δ, Z=Z)))
    end
end

E = (-8e-4:1e-6:8e-4)
lines(E, (x -> f0(x; t=0.1)).(E))
gss = (x -> g(x; Z=1, t=0.7)).(E)
lines(E, gss)

##########################################
function gs(V, E; Δ=0.0001, Z=0, μ=0, t=0, n=10)
    G0m = [G0(V, E, i; μ=0, t=0, n=n) for i in 1:4n+4]
    Qm = [Q(V, E, i, j; Δ=Δ, Z=Z, n=n) for i in 1:4n+4, j in 1:4n+4]
    Pm = [P(V, E, i; Δ=Δ, Z=Z, μ=μ, t=t, n=n) for i in 1:4n+4]
    Im = diagm(ones(4n + 4))
    prob = LinearProblem(Qm-Im, -Pm; u0=G0m)
    solve(prob, IterativeSolversJL_CG()).u
end

gs(1e-4, 1e-6; n=10)[11]

E = -4e-4:1e-6:4e-4
fs = [gs(1e-4, ϵ; Z=1, t=1)[11] for ϵ in E]
lines(E, fs)
