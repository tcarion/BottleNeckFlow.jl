function to_dict(p)
    d = Dict()
    for fn in fieldnames(typeof(p))
        push!(d, string(fn) => getfield(p, fn))
    end
    d
end

function _convert(::Type{T}, fname::String) where T
    fns = string.(fieldnames(T))
    Dataset(fname, "r") do ds
        T(Tuple(ds.attrib[fn] for fn in fns)...)
    end
end