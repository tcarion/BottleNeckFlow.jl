function beta(u::UGrid, v::VGrid, dt)
    umax = maximum(u)
    vmax = maximum(v)

    (abs(umax) + abs(vmax)) / min(u.grid.dx, u.grid.dy) * dt
end
beta(sim::Simulation) = beta(sim.u, sim.v, sim.dt)

function diagr(sim::Simulation)
    getnu(sim.params, sim.u.grid) * sim.dt / min(sim.u.grid.dx, sim.u.grid.dy)^2
end

function reynh(sim::Simulation)
    beta(sim) / diagr(sim)
end

function reynhvort(sim::Simulation)
    vort = VortGrid(sim.u, sim.v)
    maxvort = maximum(vort)
    maxvort * min(sim.u.grid.dx, sim.v.grid.dy)^2 / getnu(sim.params, sim.u.grid)
end

flowrate(v::AbstractVector, h) = trapezoid(v, h)

function flowrate(g::AbstractGrid)
    Qs = Float64[]
    for col in eachrow(truegrid(g))
        push!(Qs, flowrate(col, g.grid.dy))
    end
    Qs
end