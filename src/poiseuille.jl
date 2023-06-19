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
    for I in eachindex(grid)
        i, j = Tuple(I)
        grid[I] = poiseuille(ys[j], c, Q)
    end
end
# # takes poiseuille! x 10 time
# function poiseuille2!(grid::AbstractGrid, c::GridBox, Q::Real)
#     ys = getys(grid)
#     for i in axes(grid, 1)
#         for j in axes(grid, 2)
#             grid[i, j] = poiseuille(ys[j], c, Q)
#         end
#     end
# end

# # takes poiseuille! x 2 time
# function poiseuille3!(grid::AbstractGrid, c::GridBox, Q::Real)
#     ys = getys(grid)
#     for j in axes(grid, 2)
#         for i in axes(grid, 1)
#             grid[i, j] = poiseuille(ys[j], c, Q)
#         end
#     end
# end

# # takes poiseuille! same time
# function poiseuille4!(grid::AbstractGrid, c::GridBox, Q::Real)
#     ys = getys(grid)
#     for i in eachindex(grid)
#         grid[i] = poiseuille(ys[i[2]], c, Q)
#     end
# end

# poiseuille(y::Real, canal::CanalConfig, params::SimParam) = poiseuille(y, canal, params.Q)
poiseuille!(y::AbstractGrid, gb::GridBox, params::SimParam) = poiseuille!(y, gb, params.Q)

poiseuille!(sim::Simulation) = poiseuille!(sim.u, sim.u.grid, sim.params)
poiseuille(sim::Simulation, y) = poiseuille(y, sim.u.grid, sim.params.Q)