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

function prompeaks(x, n=1; w=1, dip=false)
    if dip
        pki = argminima(x, w)
    else
        pki = argmaxima(x, w)
    end
    if isempty(pki)
        return ()
    end
    pkp = peakproms(pki, x)[2]
    sp = sortperm(pkp, rev=true)
    m = n <= length(sp) ? n : length(sp)
    peakwidths(pki[sp[1:m]], x, pkp[sp[1:m]])
end

function findpks2(data, Ic::Val, n=2; baseline=Inf, cutp=0.6)
    m = length(data[2])
    ncl = floor(Int, cutp * m)
    ncr = floor(Int, (1 - cutp) * m)
    datal = (data[1], data[2][1:ncl], data[3][1:ncl, :])
    datar = (data[1], data[2][ncr:end], data[3][ncr:end, :])
    [findpks(datal, Ic, n ÷ 2; baseline=baseline) findpks(datar, Ic, n ÷ 2; baseline=baseline)]
end

function findpks(data, ::Val{:inner}, n=1; baseline=Inf)
    Ics = Array{Union{eltype(data[3]),Missing}}(undef, length(data[1]), n)
    for (i, v) in enumerate(data[1])
        y, z = data[2], data[3][:, i]
        pk = prompeaks(z, n)
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

function findpks(data, ::Val{:median}, n=1; baseline=Inf)
    Ics = Array{Union{eltype(data[3]),Missing}}(undef, length(data[1]), n)
    for (i, v) in enumerate(data[1])
        y, z = data[2], data[3][:, i]
        pk = prompeaks(z, n)
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

function findpks(data, ::Val{:outer}, n=1; baseline=Inf)
    Ics = Array{Union{eltype(data[3]),Missing}}(undef, length(data[1]), n)
    for (i, v) in enumerate(data[1])
        y, z = data[2], data[3][:, i]
        pk = prompeaks(z, n)
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

# function intIbias(I, R)
#     interp = LinearInterpolation(R, I)
#     V = similar(I)
#     for i in eachindex(I)
#         # V[i] = I[i] > 0 ? integral(interp, 0..I[i]) : -integral(interp, (I[i])..0)
#         V[i] = solve(IntegralProblem((x, p)->interp(x), zero(I[i]), I[i]), QuadGKJL()).u
#     end
#     V
# end

# function intIbiasmap_st(I, Rs)
#     intV = similar(Rs)
#     for i in axes(Rs, 2)
#         intV[:,i] = intIbias(I, Rs[:,i])
#     end
#     return intV
# end

# function intIbiasmap(I, Rs)
#     intV = similar(Rs)
#     cols = size(Rs, 2)
#     parts = collect(Iterators.partition(1:cols, ceil(Int, cols/Threads.nthreads())))
#     ts = map(parts) do part
#         Threads.@spawn intIbiasmap_st(I, Rs[:,part])
#     end
#     intVparts = fetch.(ts)
#     for (i, part) in enumerate(parts)
#         @views intV[:,part] = intVparts[i]
#     end
#     return intV
# end

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
    interp = LinearInterpolation(I, V; extrapolate=true)
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

sortd(x, y) = (sp = sortperm(x); (x[sp], y[sp]))
function sortd(x, y, z)
    spx = sortperm(x)
    spy = sortperm(y)
    zs = z[:, spx]
    zs = zs[spy, :]
    x[spx], y[spy], zs
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

function alignx(x, y, z; yval=y[1], order=12, w=nothing)
    yaddr = argmin(abs.(y .- yval))
    xv, yv = x, z[yaddr, :]
    isnothing(w) || (yv = hma(yv, w))
    interp = LinearInterpolation(yv, xv)
    df(x) = central_fdm(order, 1)(y -> interp(y), x)
    dfs = abs.(df.(xv))
    pki, _, _, _ = prompeaks(dfs, 2)
    p = plot(xv, dfs)
    scatter!(xv[pki], dfs[pki])
    x .- sum(x[pki]) / 2, y, z, p
end
function aligny(x, y, z; q=length(x) ÷ 4, baseline=Inf)
    Ics = loedata(x, findpks2((x, y, z), Val(:inner), baseline=baseline), q=q)
    p = heatmap(x, y, z)
    scatter!(x, Ics, ms=1)
    x, y .- sum(Ics) / length(Ics), z, p
end

function flip(x; nodes=[-Inf, Inf], rev=false)
    for i in eachindex(nodes)[1:end-1]
        if nodes[i] <= x < nodes[i+1]
            return rev ? (-1)^(i + 1) : (-1)^i
        end
    end
    return rev ? (-1)^length(nodes) : (-1)^(length(nodes) - 1)
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

function jBtojx(B, Ic, xnodesidx, xrange, W; rev=false, δB=0)
    nodes = [1; xnodesidx; length(B)]
    Icflip = Ic .* flip.(B, nodes=B[nodes]; rev=false)
    Icflip[argmin(abs.(B))] > 0 || (Icflip = Ic .* flip.(B, nodes=B[nodes]; rev=true))
    Icinterp = LinearInterpolation(Icflip, B)
    δnodes = [B[xnodesidx[i+1]] - B[xnodesidx[i]] for i in eachindex(xnodesidx)[1:end-1]]
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
    IckE = LinearInterpolation(Icinterp.(B), ks)
    gap = [0; abs.(Ic[nodes]) .* flip.(nodes, nodes=nodes; rev=rev); 0]
    IckO = LinearInterpolation(gap, ks[nodes])
    # jc(x) = (1/2π)*abs(solve(
    #     IntegralProblem((k, p) -> (IckE(k) + im*IckO(k)) * exp(-im*k*x), ks[1], ks[end]), QuadGKJL()
    #     ).u
    # )
    # jcx = jc.(xrange)
    kernel(k, x) = (IckE(k) + im * IckO(k)) * exp(-im * k * x)
    jc(x) = (1 / 2π) * abs(integrate(ks, kernel.(ks, x)))
    jcx = jc.(xrange)
    return jcx, IckE.(ks), IckO.(ks), δB
end