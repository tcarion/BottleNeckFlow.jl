function create(name::String, sim::Simulation)
    attrib = Dict{String, Any}("start_time" => string(Dates.now()))
    merge!(attrib, to_dict(sim.params))
    push!(attrib, "dt" => sim.dt)
    merge!(attrib, to_dict(sim.canal))
    merge!(attrib, to_dict(sim.u.grid))
    
    Dataset(name,"c",attrib = attrib) do ds
        defDim(ds, "steps", Inf)
        defVar(ds, "steps", [0], ("steps",))
    
        defVar(ds,"u",_3dim(sim.u),("xu","yu", "steps"), attrib = Dict(
               "units" => "m/s",
               "comments" => "horizontal velocity"
        ))
        defVar(ds, "xu", getxs(sim.u) |> collect, ("xu",))
        defVar(ds, "yu", getys(sim.u) |> collect, ("yu",))
    
        # defVar(ds,"v",sim.u,("xv","yv"), attrib = Dict(
        defVar(ds,"v",_3dim(sim.v),("xv","yv", "steps"), attrib = Dict(
            "units" => "m/s",
            "comments" => "vertical velocity"
        ))
        defVar(ds, "xv", getxs(sim.v) |> collect, ("xv",))
        defVar(ds, "yv", getys(sim.v) |> collect, ("yv",))
    
        defVar(ds,"p",_3dim(sim.p),("xc","yc", "steps"), attrib = Dict(
            "units" => "m^2/s^2",
            "comments" => "reduced pressure"
        ))
        defVar(ds, "xc", getxs(sim.p) |> collect, ("xc",))
        defVar(ds, "yc", getys(sim.p) |> collect, ("yc",))
    end
end

function add_to_nc(name::String, sim::Simulation)
    Dataset(name,"a") do ds
        writstep = length(ds["steps"]) + 1
        _add_to_var(ds, "u", sim.u, writstep)
        _add_to_var(ds, "v", sim.v, writstep)
        _add_to_var(ds, "p", sim.p, writstep)
        ds["steps"][writstep] = sim.step
    end
end

function _add_to_var(ds, var, grid, newdim)
    v = ds[var]
    v[:, :, newdim] = _3dim(grid)
end

_3dim(grid) = reshape(grid |> Matrix, (size(grid)..., 1))

function UGrid(rast::AbstractRaster)
    
end