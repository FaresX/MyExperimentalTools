function fastload(data, name; type=Float64)
    for nm in keys(data)
        if occursin(name, nm)
            return type <: Real ? replace(tryparse.(type, data[nm]), nothing => NaN) : data[nm]
        end
    end
end

function resize(z, m, n; fillms=0)
    return if length(z) > m * n
        @views reshape(z[1:m*n], m, n)
    else
        @views reshape([reshape(z, :); fill(fillms, m * n - length(z))], m, n)
    end
end

function xyzforheatmap(x, y, z)
    lx, ly = length(x), length(y)
    ny = mean(resize(y, ly ÷ lx, lx), dims=2)
    nz = resize(z, ly ÷ lx, lx)
    x[1] < x[end] || (nz = reverse(nz, dims=2))
    ny[1] < ny[end] || (nz = reverse(nz, dims=1))
    return sort(x), sort(reshape(ny, :)), nz
end

function findIcs2(data, Ic::Val, n=2; baseline=Inf, cutp=0.6)
    m = length(data[2])
    ncl = floor(Int, cutp * m)
    ncr = floor(Int, (1 - cutp) * m)
    datal = (data[1], data[2][1:ncl], data[3][1:ncl, :])
    datar = (data[1], data[2][ncr:end], data[3][ncr:end, :])
    [findIcs(datal, Ic, n ÷ 2; baseline=baseline) findIcs(datar, Ic, n ÷ 2; baseline=baseline)]
end

function findIcs(x, y, z, n=1; method=:inner, baseline=Inf, min=nothing, max=nothing)
    findIcs((x, y, z), Val(method), n; baseline=baseline, min=min, max=max)
end

function filterIcs!(x, y, z, Ic; filtersz=[], filtersy=[])
    for ft in filtersz
        for (xj, ic) in enumerate(Ic)
            yi = argmin(abs.(y .- ic))
            ft[1] < z[yi, xj] < ft[2] && (Ic[xj] = missing)
        end
    end
    for ft in filtersy
        for (xj, ic) in enumerate(Ic)
            idxl = argmin(abs.(x .- ft[1]))
            idxr = argmin(abs.(x .- ft[2]))
            x[idxl] < x[xj] < x[idxr] && !ismissing(ic) && ft[3] < ic < ft[4] && (Ic[xj] = missing)
        end
    end
    return Ic
end
filterIcs(x, y, z, Ic; filtersz=[], filtersy=[]) = filterIc!(x, y, z, copy(Ic); filtersz=filtersz, filtersy=filtersy)
function limitIcsboundsy!(f, x, Ic)
    for (i, (xi, ic)) in enumerate(zip(x, Ic))
        lb, ub = f(xi)
        (!ismissing(ic) && lb < ic < ub) || (Ic[i] = missing)
    end
    return Ic
end
limitpksboundsy(f, x, Ic) = limitIcsboundsy!(f, x, copy(Ic))
function limitIcsboundsy!(x, Ic, lbs, ubs)
    lbxs = [p[1] for p in lbs]
    lbys = [p[2] for p in lbs]
    lbxssp, lbyssp = sortd(lbxs, lbys)
    ubxs = [p[1] for p in ubs]
    ubys = [p[2] for p in ubs]
    ubxssp, ubyssp = sortd(ubxs, ubys)
    lbinterp = LinearInterpolation(lbyssp, lbxssp; extrapolate=true)
    ubinterp = LinearInterpolation(ubyssp, ubxssp; extrapolate=true)
    limitIcsboundsy!(x, Ic) do b
        lbinterp(b), ubinterp(b)
    end
    return Ic, x -> lbinterp(x), x -> ubinterp(x)
end
limitpksboundsy(x, Ic, lbs, ubs) = limitIcsboundsy!(x, copy(Ic), lbs, ubs)
function falldownIcs!(Ics)
    for i in axes(Ics, 1)
        avail = collect(eltype(Ics), skipmissing(Ics[i,:]))
        Ics[i,:] .= append!(avail, fill(missing, size(Ics, 2) - length(avail)))
    end
    return Ics
end
falldownIcs(Ics) = falldownIcs!(copy(Ics))


function findIcs(data, ::Val{:inner}, n=1; baseline=Inf, min=nothing, max=nothing)
    Ics = Array{Union{eltype(data[3]),Missing}}(undef, length(data[1]), n)
    for (i, v) in enumerate(data[1])
        y, z = data[2], data[3][:, i]
        pk = prompeaks(z, n; min=min, max=max)
        if isempty(pk)
            continue
        else
            pki, pkp, pkle, pkre = pk
        end
        dips = argminima(z)
        dipi = isempty(dips) ? argmin(z) : dips[argmin(z[dips])]
        if z[dipi] < baseline
            Icind = [y[pki[i]] < 0 ? round(Int, pkre[i]) : round(Int, pkle[i]) for i in eachindex(pki)]
            m = length(Icind)
            Ic = sort(y[Icind])
            if m < n
                Ics[i, 1:m] = Ic
            else
                Ics[i, :] = Ic
            end
        else
            Ics[i, :] = zeros(n)
        end
    end
    Ics
end

function findIcs(data, ::Val{:median}, n=1; baseline=Inf, min=nothing, max=nothing)
    Ics = Array{Union{eltype(data[3]),Missing}}(undef, length(data[1]), n)
    for (i, v) in enumerate(data[1])
        y, z = data[2], data[3][:, i]
        pk = prompeaks(z, n; min=min, max=max)
        if isempty(pk)
            continue
        else
            pki, pkp, pkle, pkre = pk
        end
        dips = argminima(z)
        dipi = isempty(dips) ? argmin(z) : dips[argmin(z[dips])]
        if z[dipi] < baseline
            Icind = pki
            m = length(Icind)
            Ic = sort(y[Icind])
            if m < n
                Ics[i, 1:m] = Ic
            else
                Ics[i, :] = Ic
            end
        else
            Ics[i, :] = zeros(n)
        end
    end
    Ics
end

function findIcs(data, ::Val{:outer}, n=1; baseline=Inf, min=nothing, max=nothing)
    Ics = Array{Union{eltype(data[3]),Missing}}(undef, length(data[1]), n)
    for (i, v) in enumerate(data[1])
        y, z = data[2], data[3][:, i]
        pk = prompeaks(z, n; min=min, max=max)
        if isempty(pk)
            continue
        else
            pki, pkp, pkle, pkre = pk
        end
        dips = argminima(z)
        dipi = isempty(dips) ? argmin(z) : dips[argmin(z[dips])]
        if z[dipi] < baseline
            Icind = [y[pki[i]] > 0 ? round(Int, pkre[i]) : round(Int, pkle[i]) for i in eachindex(pki)]
            m = length(Icind)
            Ic = sort(y[Icind])
            if m < n
                Ics[i, 1:m] = Ic
            else
                Ics[i, :] = Ic
            end
        else
            Ics[i, :] = zeros(n)
        end
    end
    Ics
end

function loedata(x, y; q=round(Int, length(x) / 4))
    if ndims(y) == 1
        return loess(x, y; q=q).(x)
    elseif ndims(y) == 2
        ny = similar(y)
        for i in 1:size(y)[2]
            ny[:, i] = loess(x, y[:, i]; q=q).(x)
        end
        return ny
    end
end

# function loedata(x, y; q=round(Int, length(data[1])/4))
#     data[1], loedata(data[1], data[2]; q=q)
# end

function loedata(x, y, z; q=round(Int, length(y) / 4), row=false)
    if row
        zt = loedata(x, transpose(z); q=q)
        return x, y, Matrix(transpose(zt))
    else
        return x, y, loedata(y, z; q=q)
    end
end

function Smoothers.hma(y::Matrix, w=5; row=false)
    yn = similar(y)
    if row
        for j in 1:size(y, 2)
            yn[:, j] = hma(y[:, j], w)
        end
    else
        for i in 1:size(y, 1)
            yn[i, :] = hma(y[i, :], w)
        end
    end
    yn
end

function transversal(x, y, z, p_value, orientation=:-)
    if orientation == :-
        p_addr = argmin(abs.(y .- p_value))
        sp_x = sortperm(x)
        return x[sp_x], z[p_addr, sp_x]
    else
        p_addr = argmin(abs.(x .- p_value))
        sp_y = sortperm(y)
        return y[sp_y], z[sp_y, p_addr]
    end
end

function intIbias(I, R::Vector)
    idx0 = argmin(abs.(I))
    cumint = cumul_integrate(I, R)
    return cumint .- cumint[idx0]
end

function intIbias(I, R::Matrix)
    idx0 = argmin(abs.(I))
    cumint = cumul_integrate(I, R)
    return cumint .- cumint[idx0, :]'
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

function determineVs(V, δ)
    minV, maxV = extrema(V)
    binV = collect(minV:δ:maxV)
    idx = argmin(abs.(binV))
    offset = binV[idx] < 0 ? (binV[idx] + binV[idx+1])/2 : (binV[idx] + binV[idx-1])/2
    binV .-= offset
    return offset < 0 ? [binV[1] - δ; binV] : [binV; binV[end] + δ]
end
binstep(V::Vector, I, δ) = binstep_kernel(V, I, determineVs(V, δ))

function binstep(Vs::AbstractMatrix, I, δ)
    binV = determineVs(Vs, δ)
    m = length(binV)
    Ihs = Matrix{Union{Missing,Float64}}(undef, m - 1, size(Vs, 2))
    for j in axes(Vs, 2)
        @views Ihs[:, j] .= binstep_kernel(Vs[:, j], I, binV)[2]
    end
    binV[2:end], Ihs
end

# function binstepmap(Vs, I, δ)
#     minV, maxV = extrema(Vs)
#     binV = collect(minV:δ:maxV)
#     m = length(binV)
#     Ihs = Matrix{Float64}(undef, m-1, size(Vs, 2))
#     for vcol in axes(Vs, 2)
#         V = Vs[:,vcol]
#         Ih = zeros(m-1)
#         minv, maxv = extrema(V)
#         for i in eachindex(binV)[1:end-1]
#             if binV[i+1] < minv || binV[i] >= maxv
#                 continue
#             elseif binV[i] <= minv && binV[i+1] > minv
#                 r = argmin(abs.(V.-binV[i+1]))
#                 Ih[i] = - -(extrema(I[1:r])...)
#             elseif binV[i] >= minv && binV[i+1] < maxv
#                 l = argmin(abs.(V.-binV[i]))
#                 r = argmin(abs.(V.-binV[i+1]))
#                 Ih[i] = - -(extrema(I[l:r])...)
#             elseif binV[i] >= minv && binV[i+1] > maxv
#                 l = argmin(abs.(V.-binV[i]))
#                 Ih[i] = - -(extrema(I[l:end])...)
#             elseif binV[i] <= minv && binV[i+1] > maxv
#                 Ih[i] = - -(extrema(I)...)
#             end
#         end
#         Ihs[:,vcol] = Ih
#     end
#     binV, Ihs
# end

function IV2dIdV(V, I, order=12)
    # interp = LinearInterpolation(I, V; extrapolate=true)
    interp = linear_interpolation(sortd(V, I)...; extrapolation_bc=Line())
    df(x) = central_fdm(order, 1)(y -> interp(y), x)
    df.(V)
end

function IV2dIdVmap(Vs, I, order=12)
    Gs = similar(Vs)
    for col in axes(Vs, 2)
        Gs[:, col] = IV2dIdV(Vs[:, col], I, order)
    end
    Gs
end

function sortd(x, y, z::Matrix)
    spx = sortperm(x)
    spy = sortperm(y)
    zs = z[:, spx]
    zs = zs[spy, :]
    x[spx], y[spy], zs
end
function sortd(x, ys...)
    sp = sortperm(x)
    return x[sp], [y[sp] for y in ys]...
end

function normalization(z)
    zn = copy(z)
    for j in axes(zn, 2)
        all(ismissing, z[:, j]) && continue
        minj, maxj = extrema(skipmissing(z[:, j]))
        maxj == minj && continue
        for i in axes(z, 1)
            zn[i, j] = (z[i, j] - minj) / (maxj - minj)
        end
    end
    zn
end

# function interpVs(Vs, Rs)
#     minv, maxv = extrema(Vs)
#     rangev = range(minv, maxv, length=size(Rs, 1))
#     Rsn = similar(Rs, Union{Float64,Missing})
#     fill!(Rsn, missing)
#     for j in axes(Rs, 2)
#         # interp = LinearInterpolation(Rs[:, j], Vs[:, j]; extrapolate=true)
#         interp = linear_interpolation(sortd(Vs[:,j], Rs[:, j])...; extrapolation_bc=Line())
#         minvj, maxvj = extrema(Vs[:, j])
#         idxmin, idxmax = argmin(abs.(rangev .- minvj)), argmin(abs.(rangev .- maxvj))
#         Rsn[idxmin:idxmax, j] = interp.(rangev[idxmin:idxmax])
#     end
#     rangev, Rsn
# end

# function alignx(x, y, z; yval=y[1], order=12, w=nothing)
#     yaddr = argmin(abs.(y.-yval))
#     xv, yv = x, z[yaddr,:]
#     isnothing(w) || (yv = hma(yv, w))
#     interp = LinearInterpolation(yv, xv)
#     df(x) = central_fdm(order,1)(y->interp(y), x)
#     dfs = abs.(df.(xv))
#     pki,_,_,_ = prompeaks(dfs, 2)
#     p = plot(xv, dfs)
#     scatter!(xv[pki], dfs[pki])
#     x.-sum(x[pki])/2, y, z, p
# end
# function aligny(x, y, z; q=length(x)÷4, baseline=Inf)
#     Ics = loedata(x, findIcs2((x, y, z), Val(:inner), baseline=baseline), q=q)
#     p = heatmap(x,y,z)
#     scatter!(x, Ics, ms=1)
#     x, y.-sum(Ics)/length(Ics), z, p
# end

function flip(x; nodes=[-Inf, Inf], rev=false)
    for i in eachindex(nodes)[1:end-1]
        if nodes[i] <= x < nodes[i+1]
            return rev ? (-1)^(i + 1) : (-1)^i
        end
    end
    n = length(nodes)
    # x == nodes[end] && return rev ? (-1)^(n+1) : (-1)^n
    return rev ? (-1)^n : (-1)^(n - 1)
end

function interpxyz(x, y, z, nrange, h=false)
    sz = copy(z)
    x[1] < x[end] || (sz = reverse(sz, dims=2))
    y[1] < y[end] || (sz = reverse(sz, dims=1))
    nx, ny = sort(x), sort(y)
    if h
        nz = Matrix{eltype(z)}(undef, size(sz, 1), length(nrange))
        for i in axes(sz, 1)
            nz[i, :] = LinearInterpolation(sz[i, :], x).(nrange)
        end
        return nz
    else
        nz = Matrix{eltype(z)}(undef, length(nrange), size(sz, 2))
        for j in axes(sz, 2)
            nz[:, j] = LinearInterpolation(sz[:, j], y).(nrange)
        end
        return nz
    end
end

function loadxy(path; xl="Ibias", yl="Vx;Ix")
    data = load(path, "data")
    x = fastload(data, xl)
    vl, il = split(yl, ';')
    y = fastload(data, vl) ./ fastload(data, il)
    return x, y
end

function loadxyz(path; xl="RF-pow", yl="Ibias", zl="Vx;Ix")
    data = load(path, "data")
    x = fastload(data, xl)
    y = fastload(data, yl)
    vl, il = split(zl, ';')
    z = fastload(data, vl) ./ fastload(data, il)
    xyzforheatmap(x, y, z)
end

# function alignIc(B, Ic, dBrange)
#     Ic .-= sum(Ic)/length(Ic)
#     s2 = []
#     for i in eachindex(dBrange)
#         Bn = B.+dBrange
#         ind0 = argmin(abs.())
# function jBtojx(x, B, Ic; nodes=nothing)
#     b = 2π*3E-7/ustrip(Unitful.Φ0)
#     Brange = B[1]..B[end]
#     abs(b*integral(B->Icsf(B)*exp(-im*b*B*x), (Bs[1]..Bs[end])*0.6E-4))
# root_path="Y:\\110-Triton2016\\DESKTOP-1EFK5BE\\D\\Exp_data"
# path="$root_path/2022-10/221028/221028.012"
# path="$root_path/2022-11/221101/221101.002"
# data = CSV.read(path, DataFrame; header=0)
# x1, x2, y, z = data[!,14]*1E3, data[!,15]*1E3, data[!,4]*0.5E9, data[!,8]./data[!,6]
# binx(x2, 0.559*0.99)
# x2[362]-x2[361]

# function jBtojx(B, Ic, xnodesidx, xrange, W; rev=false, δB=0)
#     # nodes = [1; xnodesidx; length(B)]
#     nodes = xndoexidx
#     Beff = B[nodes]
#     Iceff = Ic[nodes]
#     Icflip = Ic .* flip.(B, nodes=B[nodes]; rev=false)
#     Icflip[argmin(abs.(B))] > 0 || (Icflip = Ic .* flip.(B, nodes=B[nodes]; rev=true))
#     # Icinterp = LinearInterpolation(Icflip, B)
#     Icinterp = linear_interpolation(B, Icflip; extrapolation_bc=Line())
#     δnodes = [B[xnodesidx[i+1]] - B[xnodesidx[i]] for i in eachindex(xnodesidx)[1:end-1]]
#     if δB == 0
#         δB = if length(δnodes) == 1
#             mean(δnodes) / 2
#         else
#             maxidx = argmax(abs.(δnodes))
#             deleteat!(δnodes, maxidx)
#             mean(δnodes)
#         end
#     end
#     Λ = ustrip(Unitful.Φ0) / (δB * W)
#     ks = 2π * Λ * B / ustrip(Unitful.Φ0)
#     # IckE = LinearInterpolation(Icinterp.(B), ks)
#     IckE = linear_interpolation(ks, Icinterp.(B); extrapolation_bc=Line())
#     # gap = [0; abs.(Ic[nodes]) .* flip.(nodes, nodes=nodes; rev=rev); 0]
#     gap = abs.(Ic[nodes]) .* flip.(nodes, nodes=nodes; rev=rev)
#     # IckO = LinearInterpolation(gap, ks[nodes])
#     IckO = linear_interpolation(ks[nodes], gap; extrapolation_bc=Line())
#     # jc(x) = (1/2π)*abs(solve(
#     #     IntegralProblem((k, p) -> (IckE(k) + im*IckO(k)) * exp(-im*k*x), ks[1], ks[end]), QuadGKJL()
#     #     ).u
#     # )
#     # jcx = jc.(xrange)
#     kernel(k, x) = (IckE(k) + im * IckO(k)) * exp(-im * k * x)
#     # jc(x) = (1 / 2π) * abs(integrate(ks, kernel.(ks, x)))
#     jc(x) = (1 / 2π) * abs(solve(IntegralProblem((k, p) -> (kernel(k, x)), (B[xnodesidx[1]], B[xnodesidx[end]])), QuadGKJL()).u)
#     jcx = jc.(xrange)
#     return jcx, IckE.(ks), IckO.(ks), δB
# end

function jBtojx(B, Ic, nodes, xrange, W; rev=false, δB=0)
    Beff = B[nodes[1]:nodes[end]]
    Iceff = Ic[nodes[1]:nodes[end]]
    Icflip = Iceff .* flip.(Beff, nodes=B[nodes]; rev=false)
    Icflip[argmin(abs.(Beff))] > 0 || (Icflip = Iceff .* flip.(Beff, nodes=B[nodes]; rev=true))
    Icinterp = linear_interpolation(Beff, Icflip; extrapolation_bc=Line())
    δnodes = [B[nodes[i+1]] - B[nodes[i]] for i in eachindex(nodes)[1:end-1]]
    if δB == 0
        δB = if length(δnodes) == 1
            mean(δnodes) / 2
        else
            maxidx = argmax(abs.(δnodes))
            deleteat!(δnodes, maxidx)
            mean(δnodes)
        end
    end
    Λ = ustrip(Unitful.Φ0) / (δB * W)
    ks = 2π * Λ * B / ustrip(Unitful.Φ0)
    kseff = ks[nodes[1]:nodes[end]]
    IckE = linear_interpolation(kseff, Icinterp.(Beff); extrapolation_bc=Line())
    gap = abs.(Ic[nodes]) .* flip.(nodes, nodes=nodes; rev=rev)
    gap[end] *= -1
    IckO = linear_interpolation(ks[nodes], gap; extrapolation_bc=Line())
    kernel(k, x) = (IckE(k) + im * IckO(k)) * exp(-im * k * x)
    jc(x) = (1 / 2π) * abs(solve(IntegralProblem((k, p) -> (kernel(k, x)), kseff[1], kseff[end]), QuadGKJL()).u)
    jcx = jc.(xrange)
    return jcx, IckE.(kseff), IckO.(kseff), δB
end

function limitdata(x, y, z, limits)
    xlimits, ylimits = limits
    if isnothing(xlimits)
        xl, xr = extrema(x)
    else
        xl = isnothing(xlimits[1]) ? min(x...) : xlimits[1]
        xr = isnothing(xlimits[2]) ? max(x...) : xlimits[2]
    end
    if isnothing(ylimits)
        yd, yu = extrema(y)
    else
        yd = isnothing(ylimits[1]) ? min(y...) : ylimits[1]
        yu = isnothing(ylimits[2]) ? max(y...) : ylimits[2]
    end
    xlidx = argmin(abs.(x .- xl))
    xridx = argmin(abs.(x .- xr))
    ydidx = argmin(abs.(y .- yd))
    yuidx = argmin(abs.(y .- yu))
    return x[xlidx:xridx], y[ydidx:yuidx], z[ydidx:yuidx,xlidx:xridx]
end

function limitdata(x, y, limits)
    idxl = argmin(abs.(x .- limits[1]))
    idxr = argmin(abs.(x .- limits[2]))
    return x[idxl:idxr], y[idxl:idxr]
end

# using FiniteDifferences
# using DataInterpolations
# using Interpolations
# function binstep_kernel(V::AbstractVector, I, binV)
#     sp = sortperm(V)
#     Vsp = V[sp]
#     Isp = I[sp]
#     vi = 1
#     Ih = Vector{Union{Missing, Float64}}(undef, length(binV)-1)
#     for (i, bdv) in enumerate(binV[1:end-1])
#         Ihbin = []
#         for (j, v) in enumerate(Vsp[vi:end])
#             if bdv <= v < binV[i+1]
#                 push!(Ihbin, Isp[vi + j - 1])
#             else
#                 Ih[i] = length(Ihbin) > 1 ? - -(extrema(Ihbin)...) : missing
#                 vi = vi + j
#                 break
#             end
#         end
#         vi == length(V) && break
#     end
#     binV[2:end], Ih
# end

# binstep(V::Vector, I, δ) = binstep_kernel(V, I, collect(V[1]:δ:V[end]))

# function binstep(Vs::AbstractMatrix, I, δ)
#     minV, maxV = extrema(Vs)
#     binV = collect(minV:δ:maxV)
#     m = length(binV)
#     Ihs = Matrix{Union{Missing,Float64}}(undef, m - 1, size(Vs, 2))
#     for j in axes(Vs, 2)
#         @views Ihs[:, j] .= binstep_kernel(Vs[:, j], I, binV)[2]
#     end
#     binV[2:end], Ihs
# end

function dR(V::AbstractVector, I; p=12)
    dRm = central_fdm(p, 1)
    # interpvm = LinearInterpolation(V, I; extrapolate=true)
    interpvm = linear_interpolation(I, V; extrapolation_bc=Line())
    dRm_extra(i) = extrapolate_fdm(dRm, x -> interpvm(x), i)[1]
    return dRm_extra.(I)
end
function dR(V::AbstractMatrix, I; p=12)
    R = similar(V)
    for j in axes(V, 2)
        @views R[:,j] = dR(V[:,j], I)
    end
    return R
end

function interpVs(Vs, Rs)
    minv, maxv = extrema(Vs)
    rangev = range(minv, maxv, length=size(Rs, 1))
    Rsn = similar(Rs, Union{Float64,Missing})
    fill!(Rsn, missing)
    for j in axes(Rs, 2)
        # interp = LinearInterpolation(Rs[:, j], Vs[:, j]; extrapolate=true)
        sp = sortperm(Vs[:, j])
        interp = linear_interpolation(Vs[:, j][sp], Rs[:, j][sp]; extrapolation_bc=Line())
        minvj, maxvj = extrema(Vs[:, j])
        idxmin, idxmax = argmin(abs.(rangev .- minvj)), argmin(abs.(rangev .- maxvj))
        Rsn[idxmin:idxmax, j] = interp.(rangev[idxmin:idxmax])
    end
    rangev, Rsn
end

function interpmatrix(y, z, rate=1)
    newy = range(y[1], y[end], length=round(Int, length(y)*rate))
    newz = similar(z, length(newy), size(z, 2))
    for j in axes(z, 2)
        # interp = LinearInterpolation(z[:, j], y; extrapolate=true)
        interp = linear_interpolation(y, z[:, j]; extrapolation_bc=Line())
        @views newz[:, j] .= interp.(newy)
    end
    return newy, newz
end

function lockin!(f; α=0.01, freq=100)
    g(i, t) = f(i + α*cos(2π*freq*t))*cos(2π*freq*t)
    F(i) = 2freq*solve(IntegralProblem((t, p) -> g(i, t),0.0,1.0/freq), QuadGKJL()).u/α
end

function lockin!(is, vs; α=0.01, freq=100, m=0, T=0.1, t=Inf)
    interp = LinearInterpolation(is, vs; extrapolation_bc=Line())
    F = lockin!(x -> interp(x); α=α, freq=freq)
    return rcfilter!(F.(is); m=m, T=T, t=t)
end

function rcfilter!(ys; m=0, T=0.1, t=Inf)
    h = 1 - exp(-t/T)
    for i in 1:m
        for i in eachindex(ys)[2:end]
            ys[i] = h*ys[i] + (1-h)*ys[i-1]
        end
    end
    return ys
end