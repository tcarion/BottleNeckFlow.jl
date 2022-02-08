function convectionCenteredX(u::UGrid, vloc, i, j)
    0.5 * u[i, j] * (u[i+1, j] - u[i-1, j]) / u.grid.dx 
        + 0.5 * vloc * (u[i, j+1] - u[i, j-1]) / u.grid.dy
end

function convectionCenteredY(v::UGrid, uloc, i, j)
    0.5 * uloc * (v[i+1, j] - v[i-1, j]) / v.grid.dx
        + 0.5 * v[i, j] * (v[i, j+1] - v[i, j-1]) / v.grid.dy
end

function diffCentered(g::AbstractGridGhost, i, j)
    dx2 = 1 / (g.grid.dx * g.grid.dx)
    dy2 = 1 / (g.grid.dy * g.grid.dy)
    c = g[i, j]
    dx2 * (g[i+1, j] - 2. * c + g[i-1, j]) + dy2 * (g[i, j+1] - 2. * c + g[i, j-1])
end

function pressureX(p::AbstractGrid, i, j)
    (p[i, j] - p[i-1, j]) / p.grid.dx
end

function pressureY(p::AbstractGrid, i, j)
    (p[i, j] - p[i, j-1]) / p.grid.dy
end

function euler!(g::AbstractGridGhost, Δt, f::Function)
    for ei in eachindex(g)
        i,j = Tuple(ei)
        g[ei] = g[ei] + Δt * f(g, i, j)
    end
end

function advectionDiffusionX(g::AbstractGridGhost, i, j)
    
end