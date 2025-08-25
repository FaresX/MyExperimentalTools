function prompeaks(x, n=1; w=1, dip=false, min=nothing, max=nothing)
    if dip
        pki = argminima(x, w)
    else
        pki = argmaxima(x, w)
    end
    if isempty(pki)
        return ()
    end
    pkp = peakproms(pki, x; min=min, max=max)[2]
    sp = sortperm(pkp, rev=true)
    m = n <= length(sp) ? n : length(sp)
    peakwidths(pki[sp[1:m]], x, pkp[sp[1:m]])
end

function findpks(f, x, y, z; n=1, method=:inner, dip=false, w=1, min=nothing, max=nothing)
    pkis = fill(NaN, (n, length(x)))
    pkys = fill(NaN, (n, length(x)))
    pkzs = fill(NaN, (n, length(x)))
    for j in axes(z, 2)
        pki, pky, pkz = findpks(y, z[:, j]; method=method, n=n, dip=dip, w=w, xlims=f(x[j]), min=min, max=max)
        m = length(pki)
        pkis[1:m, j] .= pki
        pkys[1:m, j] .= pky
        pkzs[1:m, j] .= pkz
    end
    return pkis, pkys, pkzs
end
function findpks(x, y, z; n=1, method=:inner, dip=false, w=1, lbs=nothing, ubs=nothing, min=nothing, max=nothing)
    minx, maxx = extrema(x)
    miny, maxy = extrema(y)
    isnothing(lbs) && (lbs = [(minx, miny), (maxx, miny)])
    lbxs = [p[1] for p in lbs]
    lbys = [p[2] for p in lbs]
    lbxssp, lbyssp = sortd(lbxs, lbys)
    isnothing(ubs) && (ubs = [(minx, maxy), (maxx, maxy)])
    ubxs = [p[1] for p in ubs]
    ubys = [p[2] for p in ubs]
    ubxssp, ubyssp = sortd(ubxs, ubys)
    # lbinterp = LinearInterpolation(lbyssp, lbxssp; extrapolate=true)
    # ubinterp = LinearInterpolation(ubyssp, ubxssp; extrapolate=true)
    lbinterp = linear_interpolation(lbxssp, lbyssp; extrapolation_bc=Line())
    ubinterp = linear_interpolation(ubxssp, ubyssp; extrapolation_bc=Line())
    return findpks(x, y, z; n=n, method=method, dip=dip, w=w, min=min, max=max) do b
        lbinterp(b), ubinterp(b)
    end..., x -> lbinterp(x), x -> ubinterp(x)
end

function findpks(x, y; n=1, method=:inner, dip=false, w=1, xlims=(-Inf, Inf), min=nothing, max=nothing)
    findpks(x, y, Val(method), n; dip=dip, w=w, xlims=xlims, min=min, max=max)
end

function findpks(x, y, ::Val{:inner}, n=1; dip=false, w=1, xlims=(-Inf, Inf), min=nothing, max=nothing)
    pki, _, _, _ = prompeaks(y, n; dip=dip, w=w, min=min, max=max)
    pkxs = x[pki]
    pkys = y[pki]
    for i in eachindex(pki)
        xlims[1] < pkxs[i] < xlims[2] || (pkxs[i] = NaN; pkys[i] = NaN)
    end
    return pki, pkxs, pkys
end

function findpks(x, y, ::Val{:meadian}, n=1; dip=false, xlims=(-Inf, Inf), min=nothing, max=nothing)
    _, _, pkle, _ = prompeaks(y, n; dip=dip, min=min, max=max)
    return pkle, x[pkle], y[pkle]
end

function findpks(x, y, ::Val{:outter}, n=1; dip=false, xlims=(-Inf, Inf), min=nothing, max=nothing)
    _, _, _, pkre = prompeaks(y, n; dip=dip, min=min, max=max)
    return pki, x[pkre], y[pkre]
end

function falldownpks!(pks)
    for j in axes(pks, 2)
        avail = pks[:, j][findall(x -> !isnan(x), pks[:, j])]
        pks[:, j] .= append!(avail, fill(NaN, size(pks, 1) - length(avail)))
    end
    return pks
end
falldownpks(pks) = falldownpks!(copy(pks))

function falldownpks!(pks, pkz)
    for j in axes(pks, 2)
        cpkis = findall(x -> !isnan(x), pks[:, j])
        availpk = pks[:, j][cpkis]
        availz = pkz[:, j][cpkis]
        pks[:, j] .= append!(availpk, fill(NaN, size(pks, 1) - length(availpk)))
        pkz[:, j] .= append!(availz, fill(NaN, size(pks, 1) - length(availz)))
    end
    return pks, pkz
end

function filterpksz!(pks, pkz; min=-Inf, max=Inf)
    for i in eachindex(pks)
        min <= pkz[i] <= max || (pks[i] = NaN; pkz[i] = NaN)
    end
    return pks, pkz
end