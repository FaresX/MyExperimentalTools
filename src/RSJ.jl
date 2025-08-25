using Distributed
nprocs() == 1 && addprocs(7)

@everywhere using DifferentialEquations
@everywhere using DataInterpolations
@everywhere using FiniteDifferences
# using Roots
@everywhere using Integrals
@everywhere using SharedArrays
# using Unitful
# using Plots
# plotly()
using GLMakie

# f(ϕ, p, τ) = 1.1 + sin(τ) - sin(ϕ) 
# ϕtprob = ODEProblem(f, 0, (0, 100))
# ϕtsol = solve(ϕtprob, Tsit5())
# interpϕt = LinearInterpolation(ϕtsol.u, ϕtsol.t)
# ϕ(τ) = interpϕt(τ)
# # # τ0 = find_zero(τ -> ϕ(τ) - 2π, (ϕtsol.t[1], ϕtsol.t[end]))
# dϕm = central_fdm(6, 1, max_range=0.01)
# dϕm_extra(τ) = extrapolate_fdm(dϕm, ϕ, τ)[1]
# dϕs = dϕm_extra.(ϕtsol.t)
# interpdϕt = LinearInterpolation(dϕs, ϕtsol.t)
# dϕ(τ, p) = interpdϕt(τ)
# v(τ) = dϕ(τ, 0)
# vτs = v.(ϕtsol.t)
# pki = argmaxima(vτs)
# intdϕtprob = IntegralProblem(dϕ, 0, 100)
# intdϕtsol = solve(intdϕtprob, QuadGKJL())
# plot(ϕtsol)
# ϕtsol.t
# lines(ϕtsol.t, dϕs)
# scatter!(ϕtsol.t[pki], vτs[pki])
# lines(ϕtsol.t, ϕtsol.u)
@everywhere function Is(Φ, ϕ0t)
    jc(x) = 1
    ϕ(x) = ϕ0t + 2π * Φ * x
    js(ϕ) = sin(ϕ) / √(1 - 0.9sin(ϕ / 2)^2) + 0.1sin(ϕ / 2)
    prob = IntegralProblem((x, p) -> jc(x) * js(ϕ(x)), 0, 1)
    return solve(prob, QuadGKJL()).u
end

Φs = -6:0.01:6
Iϕs = [Is.(Φs, ϕ0) for ϕ0 in 0:0.1:4π]
Imaxs = [max(abs.(j)...) for j in eachrow(hcat(Iϕs...))]
lines(Φs, Imaxs)
DataInspector()

# t0 = 0
# ts = range(0, 100, length=ceil(Int, 100/1))
# ps = [1]

# for i in eachindex(ts)[2:end]
#     sw = rand() > exp(-(ts[i]-t0)/10)
#     push!(ps, sw ? -ps[i-1] : ps[i-1])
#     sw && (t0 = ts[i])
# end
# lines(ts, ps)

@everywhere function makeswseries(maxt, dt; τsw=6)
    t0 = 0
    ts = range(0, maxt, length=ceil(Int, maxt / dt))
    ps = [1]
    for i in eachindex(ts)[2:end]
        sw = rand() > exp(-(ts[i] - t0) / τsw)
        push!(ps, sw ? -ps[i-1] : ps[i-1])
        sw && (t0 = ts[i])
    end
    return LinearInterpolation(ps, ts; extrapolate=true)
end

lines(0:400, sign.(makeswseries(400, 10; τsw=10).(0:400)))

mksw = makeswseries(200, 8; τsw=10)
solve(IntegralProblem((t, p) -> sign(mksw(t)), 0, 200), QuadGKJL())

@everywhere function vmean(i, iac=0.1, ωac=0.2; Φ=0, D=0.5, αM=0.1, ϕ0=-π / 2, maxt=100, dt=1, τsw=6)
    abs(i) < 1 && (ϕ0 = asin(i))
    pτ = makeswseries(maxt, dt; τsw=τsw * exp(-2abs(i * iac)^4))
    f(ϕ, p, τ) = i + iac * sin(ωac * τ) - sin(ϕ) / √(1 - D * sin(ϕ / 2)^2) - αM * sin(ϕ / 2) * sign(pτ(τ))
    # f(ϕ, p, τ) = @. i + iac * sin(ωac * τ) - sin(ϕ)- αM * sin(ϕ / 2) * sign(sin(1τ))
    # f(ϕ, p, τ) = i + iac * sin(ωac * τ) - sin(ϕ) / √(1 - 0.6sin(ϕ / 2)^2)
    # f(ϕ, p, τ) = i + iac * sin(ωac * τ) - Is(Φ, ϕ)
    ϕtprob = ODEProblem(f, ϕ0, (0, maxt))
    ϕtsol = solve(ϕtprob, Vern9())
    (ϕtsol.u[end] - ϕtsol.u[1]) / maxt
end

is = -6:0.1:6
vs = vmean.(is, 1, 1)
lines(is, vs)
# ifft(vs)
# lines(abs.(fftshift(ifft(vs))))
# @everywhere function vmean(i, iac=0.1, ωac=0.2; αM=0.1, ϕ0=-π / 2, maxt=100)
#     abs(i) < 1 && (ϕ0 = asin(i))
#     # f(ϕ, p, τ) = i + iac * sin(ωac * τ) - sin(ϕ) / √(1 - 0.6sin(ϕ / 2)^2) - αM * sin(ϕ / 2)
#     f(ϕ, p, τ) = i + iac * sin(ωac * τ) - sin(ϕ)
#     # f(ϕ, p, τ) = i + iac * sin(ωac * τ) - sin(ϕ) / √(1 - 0.6sin(ϕ / 2)^2)
#     ϕtprob = ODEProblem(f, ϕ0, (0, maxt))
#     ϕtsol = solve(ϕtprob, Tsit5())
#     interpϕt = LinearInterpolation(ϕtsol.u, ϕtsol.t)
#     ϕ(τ) = interpϕt(τ)
#     dϕm = central_fdm(6, 1, max_range=0.01)
#     dϕm_extra(τ) = extrapolate_fdm(dϕm, ϕ, τ)[1]
#     dϕs = dϕm_extra.(ϕtsol.t)
#     # interpdϕt = LinearInterpolation(dϕs, ϕtsol.t)
#     # vτs = interpdϕt.(ϕtsol.t)
#     pki = argmaxima(dϕs)
#     return if length(pki) > 4
#         (dϕs[pki[end]]-dϕs[pki[1]])/(ϕtsol.t[pki[end]]-ϕtsol.t[pki[1]])
#     else
#         (ϕtsol.u[end] - ϕtsol.u[1]) / maxt
#     end
# end

@everywhere function dR(is, iac=1, ωac=0.2; Φ=0, D=0.5, αM=0.1, maxt=100, dt=1, τsw=6)
    dRm = central_fdm(6, 1)
    # vmeans = vmean.(is, iac, ωac; Φ=Φ, D=D, αM=αM, maxt=maxt, dt=dt, τsw=τsw)
    vmeans = [vmean(is[i], iac, ωac; Φ=Φ, D=D, αM=αM, maxt=maxt, dt=dt, τsw=τsw[i]) for i in eachindex(is)]
    interpvm = LinearInterpolation(vmeans, is; extrapolate=true)
    dRm_extra(i) = extrapolate_fdm(dRm, x -> interpvm(x), i)[1]
    return vmeans, dRm_extra.(is)
end

function interpVs(Vs, Rs)
    minv, maxv = extrema(Vs)
    rangev = range(minv, maxv, length=size(Rs, 1))
    Rsn = similar(Rs, Union{Float64,Missing})
    fill!(Rsn, missing)
    for j in axes(Rs, 2)
        interp = LinearInterpolation(Rs[:, j], Vs[:, j]; extrapolate=true)
        minvj, maxvj = extrema(Vs[:, j])
        idxmin, idxmax = argmin(abs.(rangev .- minvj)), argmin(abs.(rangev .- maxvj))
        Rsn[idxmin:idxmax, j] = interp.(rangev[idxmin:idxmax])
    end
    rangev, Rsn
end

function binstep_kernel(V::AbstractVector, I, binV)
    δV = binV[2] - binV[1]
    binVidx = []
    for i in eachindex(binV)
        d, idx = findmin(abs.(V .- binV[i]))
        push!(binVidx, d < δV ? idx : missing)
    end
    Ih = [
        if (ismissing(binVidx[i]) || ismissing(binVidx[i+1]))
            missing
        else
            - -(extrema(I[binVidx[i]:binVidx[i+1]]; init=(I[binVidx[i]], I[binVidx[i+1]]))...)
        end
        for i in eachindex(binVidx)[1:end-1]
    ]
    binV[2:end], Ih
end

binstep(V::Vector, I, δ) = binstep_kernel(V, I, collect(V[1]:δ:V[end]))


function binstep(Vs::AbstractMatrix, I, δ)
    minV, maxV = extrema(Vs)
    binV = collect(minV:δ:maxV)
    m = length(binV)
    Ihs = Matrix{Union{Missing,Float64}}(undef, m - 1, size(Vs, 2))
    for j in axes(Vs, 2)
        @views Ihs[:, j] .= binstep_kernel(Vs[:, j], I, binV)[2]
    end
    binV[2:end], Ihs
end

# is = -6:0.01:6
# @time vmeans = vmean.(is, 0.56; αM=0, maxt=1000)
# @time dRs = dR(is, 0.56; αM=0, maxt=100)
# plot(is, vmeans)
# plot(is, dRs)
is = -1:0.01:1
iacs = 0.0:0.04:4
ωac = 0.2
Rs = Matrix{Float64}(undef, length(is), length(iacs))
Rss = SharedArray(Rs)
Vs = Matrix{Float64}(undef, length(is), length(iacs))
Vss = SharedArray(Vs)
t = @distributed for i in eachindex(iacs)
    v, r = dR(is, iacs[i], ωac; Φ=1, D=0.6, αM=0.3, maxt=200, dt=1, τsw=[1e12exp(-10abs(i)) for i in is])
    @views Vss[:, i] = v
    @views Rss[:, i] = r
end
@time wait(t)
fig = Figure(fontsize=36)
ax = Axis(
    fig[1, 1],
    # title=rich(rich("D", font=:italic), rich("=0.6 ω", font=:regular), subscript(rich("ac", font=:regular)), rich("=$ωac", font=:regular)),
    xlabel=rich(rich("i", font=:italic), subscript("ac")),
    ylabel=rich("i", font=:italic),
)
p = heatmap!(ax, iacs, is, Rss', colorrange=(0, max(Rss...)))
Colorbar(fig[1, 2], p, label=rich("d", rich("v", font=:italic), "/d", rich("i", font=:italic)))
fig

# V, R = interpVs(Vss, Rss)
# fig = Figure(fontsize=24)
# ax = Axis(
#     fig[1, 1],
#     title=rich(rich("D", font=:italic), rich("=0.6 ω", font=:regular), subscript(rich("ac", font=:regular)), rich("=$ωac", font=:regular)),
#     xlabel=rich(rich("i", font=:italic), subscript("ac")),
#     ylabel=rich(rich("v", font=:italic), "(ω", subscript("ac"), ")")
# )
# p = heatmap!(ax, iacs, V / ωac, R', colorrange=(0, max(Rss...)), colorscale=Makie.pseudolog10)
# Colorbar(fig[1, 2], p, label=rich("d", rich("v", font=:italic), "/d", rich("i", font=:italic)))
# DataInspector()
# fig
Vh, Ihs = binstep(Vss./ωac, is, 0.2)
fig = Figure(fontsize=24)
ax = Axis(
    fig[1, 1],
    title=rich(rich("D", font=:italic), rich("=0.6 ω", font=:regular), subscript(rich("ac", font=:regular)), rich("=$ωac", font=:regular)),
    xlabel=rich(rich("i", font=:italic), subscript("ac")),
    ylabel=rich(rich("v", font=:italic), "(ω", subscript("ac"), ")")
)
p = heatmap!(ax, iacs, Vh, Ihs', colorrange=(0, max(Rss...)), colorscale=Makie.pseudolog10)
Colorbar(fig[1, 2], p, label=rich("d", rich("v", font=:italic), "/d", rich("i", font=:italic)))
fig



### B ###
is = -4:0.02:4
iac = 1.5
ωac = 0.4
Φs = 0:0.1:4
Rs = Matrix{Float64}(undef, length(is), length(Φs))
Rss = SharedArray(Rs)
Vs = Matrix{Float64}(undef, length(is), length(Φs))
Vss = SharedArray(Vs)
t = @distributed for i in eachindex(Φs)
    v, r = dR(is, iacs[i], ωac; Φ=Φs[i], αM=0.2, maxt=400)
    @views Vss[:, i] = v
    @views Rss[:, i] = r
end
@time wait(t)
fig = Figure(fontsize=36)
ax = Axis(
    fig[1, 1],
    # title=rich(rich("D", font=:italic), rich("=0.6 ω", font=:regular), subscript(rich("ac", font=:regular)), rich("=$ωac", font=:regular)),
    xlabel=rich("Φ", font=:italic),
    ylabel=rich("i", font=:italic),
)
p = heatmap!(ax, Φs, is, Rss', colorrange=(0, max(Rss...)))
Colorbar(fig[1, 2], p, label=rich("d", rich("v", font=:italic), "/d", rich("i", font=:italic)))
fig

V, R = interpVs(Vss, Rss)
fig = Figure(fontsize=36)
ax = Axis(
    fig[1, 1],
    # title=rich(rich("D", font=:italic), rich("=0.6 ω", font=:regular), subscript(rich("ac", font=:regular)), rich("=$ωac", font=:regular)),
    xlabel=rich("Φ", font=:italic),
    ylabel=rich(rich("v", font=:italic), "(ω", subscript("ac"), ")")
)
p = heatmap!(ax, Φs, V / ωac, R', colorrange=(0, max(Rss...)), colorscale=Makie.pseudolog10)
Colorbar(fig[1, 2], p, label=rich("d", rich("v", font=:italic), "/d", rich("i", font=:italic)))
DataInspector()
fig