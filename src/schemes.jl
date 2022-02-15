function convectionCenteredX(u::AbstractGridGhost, vloc, i, j)
    0.5 * u[i, j] * (u[i+1, j] - u[i-1, j]) / u.grid.dx 
        + 0.5 * vloc * (u[i, j+1] - u[i, j-1]) / u.grid.dy
end

function convectionCenteredY(v::AbstractGridGhost, uloc, i, j)
    0.5 * uloc * (v[i+1, j] - v[i-1, j]) / v.grid.dx
        + 0.5 * v[i, j] * (v[i, j+1] - v[i, j-1]) / v.grid.dy
end

function laplacianCentered(g::AbstractGrid, i, j)
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

function fillghost!(g::UGrid, f1, f2)
    ghost1 = g.ghost[1]
    ghost2 = g.ghost[2]
    for i in eachindex(ghost1)
        ghost1[i] = f1(g[i, 1], g[i, 2], g[i, 3])
        ghost2[i] = f2(g[i, end], g[i, end-1], g[i, end-2])
    end
end

function fillghost!(g::VGrid, f1, f2)
    ghost1 = g.ghost[1]
    ghost2 = g.ghost[2]
    for i in eachindex(ghost1)
        ghost1[i] = f1(g[1, i], g[2, i], g[1, i])
        ghost2[i] = f2(g[end, i], g[end-1, i], g[end-2, i])
    end
end

function noslip(arg::Vararg{Real})
    u1, u2, u3 = arg
    -0.2 * (u3 - 5. * u2 + 15. * u1)
end

function naturalbound(arg::Vararg{Real})
    u1 = arg[1]
    u1
end

function veldiv(u::UGrid, v::VGrid, i, j)
    (u[i+1, j] - u[i, j]) / u.grid.dx + (v[i, j+1] - v[i, j]) / v.grid.dy
end

function step_AB2!(sim::Simulation, u_n_1::UGrid, v_n_1::VGrid)
    u_n = copy(sim.u)
    v_n = copy(sim.v)
    nu = getnu(sim.params, sim.u.grid)
    xu = getxs(sim.u)
    yu = getys(sim.u)

    xv = getxs(sim.v)
    yv = getys(sim.v)
    for i in 2:size(sim.u)[1]-1
        for j in 1:size(sim.u)[2]
    # for ci in eachindex(inner_u)
            # (i, j) = Tuple(ci)
            vloc = projectionU(v_n, i, j)
            vloc_n_1 = projectionU(v_n_1, i, j)
            # f1 = - 0.5 * (3. * convectionCenteredX(u_n, vloc, i, j) - convectionCenteredX(u_n_1, vloc_n_1, i, j))
            f1 = - AB2(convectionCenteredX(u_n, vloc, i, j), convectionCenteredX(u_n_1, vloc_n_1, i, j))
            f2 = - pressureX(sim.p, i, j)
            f3 = nu*laplacianCentered(u_n, i, j)

            with_logger(sim.logger) do
                @debug "f1: $f1 - f2: $f2 - f3: $f3"
            end

            mask = inbump(sim.canal, xu[i], yu[j])
            # sim.u[ci] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
            fn1 = u_n[i, j] + sim.dt * (f1 + f2 + f3)
            if mask
                sim.u[i, j] = fn1 / (1 + sim.params.dtdτ * ramp_khi(sim))
            else
                sim.u[i, j] = fn1
            end
            # sim.u[i, j] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
        end
    end

    # for ci in eachindex(inner_v)
    for i in 1:size(sim.v)[1]
        for j in 2:size(sim.v)[2]-1
        # (i, j) = Tuple(ci)
            uloc = projectionV(u_n, i, j)
            uloc_n_1 = projectionV(u_n_1, i, j)
            f1 = - AB2(convectionCenteredY(v_n, uloc, i, j), convectionCenteredY(v_n_1, uloc_n_1, i, j))
            f2 = - pressureY(sim.p, i, j)
            f3 = nu*laplacianCentered(v_n, i, j)

            mask = inbump(sim.canal, xv[i], yv[j])

            fn1 = v_n[i, j] + sim.dt * (f1 + f2 + f3)
            if mask
                sim.v[i, j] = fn1 / (1 + sim.params.dtdτ * ramp_khi(sim))
            else
                sim.v[i, j] = fn1
            end
            # sim.v[ci] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
            # sim.v[i, j] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
        end
    end

    boundaries!(sim)
end

function step_euler!(sim::Simulation)
    u_n = copy(sim.u)
    v_n = copy(sim.v)
    nu = getnu(sim.params, sim.u.grid)
    xu = getxs(sim.u)
    yu = getys(sim.u)

    xv = getxs(sim.v)
    yv = getys(sim.v)
    for i in 2:size(sim.u)[1]-1
        for j in 1:size(sim.u)[2]
    # for ci in eachindex(inner_u)
            # (i, j) = Tuple(ci)
            vloc = projectionU(v_n, i, j)
            f1 = - convectionCenteredX(u_n, vloc, i, j)
            f2 = - pressureX(sim.p, i, j)
            f3 = nu*laplacianCentered(u_n, i, j)

            with_logger(sim.logger) do
                @debug "f1: $f1 - f2: $f2 - f3: $f3"
            end

            mask = inbump(sim.canal, xu[i], yu[j])
            # sim.u[ci] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
            fn1 = u_n[i, j] + sim.dt * (f1 + f2 + f3)
            if mask
                sim.u[i, j] = fn1 / (1 + sim.params.dtdτ * ramp_khi(sim))
            else
                sim.u[i, j] = fn1
            end
            # sim.u[i, j] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
        end
    end

    # for ci in eachindex(inner_v)
    for i in 1:size(sim.v)[1]
        for j in 2:size(sim.v)[2]-1
        # (i, j) = Tuple(ci)
            uloc = projectionV(u_n, i, j)
            f1 = - convectionCenteredY(v_n, uloc, i, j)
            f2 = - pressureY(sim.p, i, j)
            f3 = nu*laplacianCentered(v_n, i, j)

            mask = inbump(sim.canal, xv[i], yv[j])

            fn1 = v_n[i, j] + sim.dt * (f1 + f2 + f3)
            if mask
                sim.v[i, j] = fn1 / (1 + sim.params.dtdτ * ramp_khi(sim))
            else
                sim.v[i, j] = fn1
            end
            # sim.v[ci] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
            # sim.v[i, j] = mask ? (f1 + f2 + f3) / (1 + sim.params.dtdτ) : (f1 + f2 + f3)
        end
    end

    boundaries!(sim)
end

function rightboundary!(sim::Simulation)
    yu = getys(sim.u)
    for col in eachcol(sim.u)
        j = col.indices[2]
        col[end] = col[end] - sim.dt * poiseuille(sim, yu[j]) * (col[end] - col[end-1]) / sim.u.grid.dx
    end
end

function boundaries!(sim::Simulation)
    fillghost!(sim.u, noslip, noslip)
    fillghost!(sim.v, noslip, naturalbound)
    rightboundary!(sim)
end

# TO BE VERIFIED, SEEMS NOT GOOD
function pgauss(p::PGrid, i, j, dx2, dy2)
    right = i != p.grid.m ? p[i+1, j] : -p[i, j]
    left = i != 1 ? p[i-1, j] : p[i, j]
    top = j != p.grid.n ? p[i, j+1] : p[i, j] # DP is 0 at walls
    bot = j != 1 ? p[i, j-1] : p[i, j] # DP is 0 at walls

    (right + left) / dx2 + (top + bot) / dy2
end

function gauss_seidel!(p::PGrid, u::UGrid, v::VGrid, α, dt)
    dx2 = u.grid.dx * u.grid.dx
    dy2 = u.grid.dy * u.grid.dy
    coef = 0.5 * dx2 * dy2 / (dx2 + dy2)
    # innerp = @view p[2:end-1, 2:end-1]
    for ci in eachindex(p)
        i, j = Tuple(ci)
        div = - veldiv(u, v, i, j) / dt
        phit = pgauss(p, i, j, dx2, dy2)
        phistar = coef * ( div + phit )
        p[ci] = α * phistar + (1 - α) * p[ci]
    end
end
gauss_seidel!(sim::Simulation) = gauss_seidel!(sim.p, sim.u, sim.v, sim.params.alpha, sim.dt)

function residual(p::PGrid, u::UGrid, v::VGrid, dt)
    R = 0
    # SHOULD TAKE INTO ACCOUNT BOUNDARIES
    is, js = inrange(p)
    for i in is
        for j in js
            # i, j = Tuple(ci)
            lap = laplacianCentered(p, i, j)
            div = veldiv(u, v, i, j)
            R += (lap - div / dt)^2
        end
    end
    R * p.grid.dx * p.grid.dy
end
function residual(sim::Simulation)
    Rs = residual(sim.p, sim.u, sim.v, sim.dt)
    grid = sim.u.grid
    sim.dt * grid.D / getum(sim.params, grid) * sqrt( 1 / (grid.L * grid.D) * Rs )
end

function step_poisson!(sim::Simulation, tol = 1e-3)
    e = 1000.
    n = 0
    maxn = 1000
    while (abs(e) > tol) && (n < maxn)
        gauss_seidel!(sim)
        e = residual(sim)
        n += 1
    end

    if n == maxn
        with_logger(sim.logger) do
            @warn "Gauss Seidel exceeded $maxn iterations! Aborted"
        end
    end

    with_logger(sim.logger) do
        @info "Poisson equation solved after $n iterations"
    end
end

function poisson_project!(u::UGrid, v::VGrid, p::PGrid, dt)
    ustar = copy(u)
    vstar = copy(v)

    is, js = inrange(u)
    for i in is
        for j in js
            u[i, j] = ustar[i, j] - dt / u.grid.dx * (p[i-1, j] - p[i, j])
        end
    end

    is, js = inrange(v)
    for i in is
        for j in js
            v[i, j] = vstar[i, j] - dt / u.grid.dy * (p[i, j] - p[i, j-1])
        end
    end
end
poisson_project!(sim::Simulation) = poisson_project!(sim.u, sim.v, sim.p, sim.dt)
# function euler!(g::AbstractGridGhost, Δt, f::Function)
#     for ei in eachindex(g)
#         i,j = Tuple(ei)
#         g[ei] = g[ei] + Δt * f(g, i, j)
#     end
# end

function euler(u_n, Δt, f_n)
    u_n + Δt * f_n
end

# function AB2(u_n, Δt, f_n, f_n_1)
#     u_n + 0.5 * Δt * (-f_n_1 + 3. * f_n)
# end

function AB2(f_n, f_n_1)
    0.5 * (-f_n_1 + 3. * f_n)
end

function advectionDiffusionX(g::AbstractGridGhost, i, j)
    
end

function ramp_khi(ramp::Int, n::Int)
    coef = n / ramp
    coef > 1. ? 1. : coef
end

ramp_khi(sim::Simulation) = ramp_khi(sim.params.ramp, sim.step)