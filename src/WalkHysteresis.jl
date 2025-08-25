mutable struct WalkHysteresis
    x::Float64
    ysw::Float64
    yre::Float64
    fsw::Function
    fre::Function
    state::Bool
    tol::Float64
end

function (w::WalkHysteresis)(dx)
    w.x += dx
    ysw = w.fsw(w.x)
    yre = w.fre(w.x)
    if abs(w.ysw - w.yre) < w.tol
        w.state = dx > 0
    else
        abs(ysw - yre) < w.tol && (w.state = dx > 0)
    end
    w.ysw = ysw
    w.yre = yre
    return w.state ? ysw : yre
end

function lockinhys!(i, state, fsw, fre; α=0.01, freq=100, tol=0.01, lt=100)
    ts = 0:1/lt/freq:1/freq
    i0 = i + α
    walk = WalkHysteresis(i0, fsw(i0), fre(i0), fsw, fre, state, tol)
    ws = [walk.state ? walk.ysw : walk.yre]
    for i in eachindex(ts)[2:end]
        push!(ws, walk(α * (cos(2π * freq * ts[i]) - cos(2π * freq * ts[i-1]))))
    end
    return 2freq * integrate(ts, ws .* cos.(2π * freq * ts)) / α
end

function lockinhys!(is, vssw, vsre; α=0.01, freq=100, m=0, T=0.1, t=Inf, initstate=true, tol=0.01, lt=100)
    interpsw = LinearInterpolation(is, vssw; extrapolation_bc=Line())
    interpre = LinearInterpolation(is, vsre; extrapolation_bc=Line())
    fsw(i) = interpsw(i)
    fre(i) = interpre(i)
    rs = map(i -> lockinhys!(i, initstate, fsw, fre; α=α, freq=freq, tol=tol, lt=lt), is)
    return rcfilter!(rs; m=m, T=T, t=t)
end