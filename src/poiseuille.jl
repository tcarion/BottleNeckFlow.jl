# function poiseuille(y::Real, canal::CanalConfig, μ, dpdx)
#     h = canal.H
#     vmax = -1 / μ * dpdx * h^2 / 2
#     vmax * (1 - (y / h)^2)
# end

# poiseuille(y::AbstractVector, c::CanalConfig, μ, dpdx) = poiseuille.(y, Ref(c), Ref(μ), Ref(dpdx))

function poiseuille(y::Real, gb::GridBox, Q)
    uc = 1.5 * Q / gb.D
    uc * (1 - (y / (0.5*gb.D))^2)
end

poiseuille(y::AbstractVector, c::GridBox, Q) = poiseuille.(y, Ref(c), Ref(Q))

function poiseuille!(grid::AbstractGrid, c::GridBox, Q::Real)
    ys = getys(grid)
    for i in eachindex(grid)
        grid[i] = poiseuille(ys[i[2]], c, Q)
    end
end

# poiseuille(y::Real, canal::CanalConfig, params::SimParam) = poiseuille(y, canal, params.Q)
poiseuille!(y::AbstractGrid, gb::GridBox, params::SimParam) = poiseuille!(y, gb, params.Q)

poiseuille!(sim::Simulation) = poiseuille!(sim.u, sim.u.grid, sim.params)
poiseuille(sim::Simulation, y) = poiseuille(y, sim.u.grid, sim.params.Q)