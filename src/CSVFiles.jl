function header(path; t=100)
    file = open(path)
    s = readline(file)
    c = 1
    while s != "" && c < t
        s = readline(file)
        c += 1
    end
    close(file)
    c
end

function fastloaddf(data, name; type=Float64)
    for nm in names(data)
        if occursin(name, nm)
            if eltype(data[!, nm]) <: AbstractString
                return type <: Real ? replace(tryparse.(type, data[!, nm]), nothing => NaN) : data[!, nm]
            else
                return data[!, nm]
            end
        end
    end
end