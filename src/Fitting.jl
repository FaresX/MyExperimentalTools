function polyfitting(xs, ys; order=1)
    f(x; p=ones(order + 1)) = sum(p[i] * x^(i - 1) for i in 1:order+1)
    @. model(x, p) = f(x; p=p)
    sp = sortperm(xs)
    fitting = curve_fit(model, xs[sp], ys[sp], ones(order + 1))
    p = fitting.param
    return model, p
end

function nonlinearfitting(model, xs, ys, p0; lower=nothing, upper=nothing)
    fmodel(p) = sum(abs2.(model(xs, p) .-  ys))
    return if !isnothing(lower) && !isnothing(upper)
        optimize(fmodel, lower, upper, p0)
    else
        optimize(fmodel, p0)
    end
end