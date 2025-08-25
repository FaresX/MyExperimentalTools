function symmetrizeRB(B, R)
    Bsp, Rsp = sortd(B, R)
    Rinterp = LinearInterpolation(Rsp, Bsp; extrapolate=true)
    Rxx(B) = (Rinterp(B) + Rinterp(-B)) / 2
    Rxy(B) = (Rinterp(B) - Rinterp(-B)) / 2
    return Bsp, Rxx.(Bsp), Rxy.(Bsp)
end
SdHdecomposition(B, R; method=:poly, kwargs...) = SdHdecomposition(B, R, Val(method); kwargs...)
function SdHdecomposition(B, R, ::Val{:poly}; order=3)
    model, p = polyfitting(B, R; order=order)
    Bsp, Rsp = sortd(B, R)
    fitR = model(Bsp, p)
    SdHR = Rsp - fitR
    return Bsp, fitR, SdHR, p
end

function SdHdecomposition(B, R, ::Val{:loess}; q=round(Int, length(B) / 4))
    Bsp, Rsp = sortd(B, R)
    p = loess(Bsp, Rsp; q=q)
    fitR = p.(Bsp)
    SdHR = Rsp - fitR
    return Bsp, fitR, SdHR, p
end

function splitSdH(invB, R)
    invBsp, Rsp = sortd(invB, R)
    idx = findfirst(x -> x > 0, invBsp)
    invBl = invBsp[1:idx-1]
    invBr = invBsp[idx:end]
    Rspl = Rsp[1:idx-1]
    Rspr = Rsp[idx:end]
    invBl, invBr, Rspl, Rspr
end

function analysisf(x, y; fs=(length(x)-1)/(- -(extrema(x)...)), norm=true, half=true, shift=true, addzero=0)
    xmin, xmax = extrema(x)
    xn, yn = preinterpdata(x, y; fs=fs)
    if addzero != 0
        zn = round(Int, addzero/2)
        xn = xmin - zn/fs:1/fs:xmax+zn/fs
        yn = [zeros(zn); yn; zeros(zn)]
    end
    n = length(xn)
    yfft = fft(yn)
    norm && (yfft ./= n/2; yfft[1] /= 2)
    fxs = fftfreq(n, fs)
    if half
        m = Int(iseven(n) ? n/2 : (n-1)/2)
        return fxs[1:m], shift ? fftshift(yfft)[end-m:end-1] : yfft[end-m:end-1]
    else
        if shift
            m = Int(iseven(n) ? n/2 : (n-1)/2)
            fxs .-= m*fs/n
            return fxs, fftshift(yfft)
        else
            return fxs, yfft
        end
    end
end

function preinterpdata(x, y; fs=(length(x)-1)/(- -(extrema(x)...)))
    xsp, ysp = sortd(x, y)
    interp = LinearInterpolation(ysp, xsp; extrapolate=true)
    xmin, xmax = extrema(x)
    xn = xmin:1/fs:xmax
    yn = interp.(xn)
    return collect(xn), yn
end

function window(f, x, y, args...; fs=(length(x)-1)/(- -(extrema(x)...)), padding=0, kwargs...)
    xn, yn = preinterpdata(x, y; fs=fs)
    append!(xn, xn[end]+1/fs:1/fs:xn[end]+padding/fs)
    append!(yn, zeros(padding))
    return xn, yn .* f(length(yn), args...; padding=padding, kwargs...)
end