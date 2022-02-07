function poiseuille(y::Real, canal::CanalConfig, μ, dpdx)
    h = canal.H
    vmax = -1 / μ * dpdx * h^2 / 8
    vmax * (y / h - (y / h)^2)
end

poiseuille(y::AbstractVector, c::CanalConfig, μ, dpdx) = poiseuille.(y, Ref(c), Ref(μ), Ref(dpdx))

function poiseuille!(grid::AbstractGrid, c::CanalConfig, μ, dpdx)
    ys = getys(grid)
    for i in eachindex(grid)
        grid[i] = poiseuille(ys[i[2]], c, μ, dpdx)
    end
end