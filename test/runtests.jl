using BottleNeckFlow
using CairoMakie
using GLMakie
using Test
using Logging
using Rasters
using DimensionalData
bn = BottleNeckFlow

ENV["JULIA_DEBUG"] = BottleNeckFlow
ENV["JULIA_DEBUG"] = ""

H = 1.
canal = CanalConfig(H)
params = SimParam()
params.cfl = 0.3
params.tol = 1e-3
params.dtdτ = 0.1
params.ramp = 400
gridbox = GridBox(canal, .1)
gridbox = GridBox(canal, .05)
Re = 250
ubar = 3.
nu = ubar * canal.L / Re
rho = 1000
mu = nu * rho

canobs = Observable(CanalConfig(2.))
x = range(0, canal.L, length = 101)
xobs = range(0, L, length = 50) |> Observable
y = δ(x, canal)

u = UGrid(gridbox)
xu = getxs(u) |> collect
yu = getys(u) |> collect
dimarray_u = DimArray(u)

v = VGrid(gridbox)
getxs(v) |> collect
getys(v) |> collect

p = PGrid(gridbox)
getxs(p) |> collect
getys(p) |> collect

poiseuille!(u, canal, params)
poiseuille!(u, gridbox, params)
# @time poiseuille!(u, gridbox, params)
# @btime poiseuille!(u, gridbox, 1)
# @btime bn.poiseuille2!(u, gridbox, 1)
# @btime bn.poiseuille3!(u, gridbox, 1)
# @btime bn.poiseuille4!(u, gridbox, 1)

vort = VortGrid(u.grid)
vort = VortGrid(u, v)
vort = VortGrid(sim.u, sim.v)

io = open("logcfl03.log", "w")

sim = Simulation(params, gridbox, canal, SimpleLogger(io))
sim = Simulation(params, gridbox, canal)

poiseuille!(sim)

# fillghost!(sim.u, noslip, noslip)
fillghost!(sim.u)
# fillghost!(sim.v, noslip, naturalbound)
fillghost!(sim.v)
bn.rightboundary!(sim)
step_euler!(sim)
for i in 1:1000 step_euler!(sim); sim.step = sim.step + 1; end
step_poisson!(sim, 1e-5)
poisson_project!(sim)

bn.boundaries!(sim)
bn.make_step!(sim)

runsim!(sim, 50000)
sim = @async runsim(gridbox, canal, 500)

xu = getxs(sim.u) |> collect
yu = getys(sim.u) |> collect
pos = [(x, y) for x in xu for y in yu]

step_euler!(sim); simobs[] = sim

gauss_seidel!(sim)
residual(sim)

simobs = Observable(Simulation(params, gridbox, canal))
uobs = Observable(simobs[].u)

Makie.Observables.onany(simobs) do sim
    uobs[] = sim.u
end
# poiseuille!(simobs[])
# step_euler!(simobs[])

divergence = DivGrid(sim.u, sim.v)

# f = Figure()
# ax = f[1, 1] = Axis(f, aspect = DataAspect())
# ax2 = f[2, 1] = Axis(f, aspect = DataAspect())

# xlims!(ax, (0, gridbox.L))
# ylims!(ax, (-0.5*gridbox.D, 0.5*gridbox.D))
# xlims!(ax2, (0, gridbox.L))
# ylims!(ax2, (-0.5*gridbox.D, 0.5*gridbox.D))

for ax in axes
    canalplot!(ax, x, canal)
end
canalplot!(ax, xobs, canobs)
canalplot(x, canal)

canobs[] = CanalConfig(1.5)
xobs[] = range(0, L-2, length= 20)

contourf!(ax, xu, yu, u, colormap = :jet)
contourf(xu, yu, u, colormap = :jet)
cont = contourf!(ax, getxs(vort), getys(vort), vort .|> abs, colormap = :jet)
cont = contourf(getxs(vort), getys(vort), vort .|> abs, colormap = :jet)
Colorbar(f[1, 2], cont)


plotgrid!(ax, sim.v, marker = :rect)
plotgrid!(ax, sim.u, marker = :circle, color = :green)

function newfig(gridbox)
    f = Figure()
    axes = [Axis(f[i, 1], aspect = DataAspect()) for i in 1:3]
    for ax in axes
        xlims!(ax, (0, gridbox.L))
        ylims!(ax, (-0.5*gridbox.D, 0.5*gridbox.D))
    end
    # ax1 = f[1, 1] = Axis(f, aspect = DataAspect(), xlims = (0, gridbox.L), ylims = (-0.5*gridbox.D, 0.5*gridbox.D))
    # ax2 = f[2, 1] = Axis(f, aspect = DataAspect(), xlims = (0, gridbox.L), ylims = (-0.5*gridbox.D, 0.5*gridbox.D))
    f, axes
end

function plotsim(axes, sim::Simulation)
    ax1, ax2, ax3 = axes
    griddedplot!(ax1, sim.u)
    griddedplot!(ax2, sim.v)
    griddedplot!(ax3, sim.p)
end
function plotsim(axes, g1::AbstractGrid, g2::AbstractGrid, g3::AbstractGrid)
    ax1, ax2, ax3 = axes
    griddedplot!(ax1, g1)
    griddedplot!(ax2, g2)
    griddedplot!(ax3, g3)
end
function plotsim(fig::Figure, g1::AbstractGrid, g2::AbstractGrid, g3::AbstractGrid)
    ax1, ax2, ax3 = axes
    griddedplot!(ax1, g1)
    griddedplot!(ax2, g2)
    griddedplot!(ax3, g3)
end
function plotsim(sim::Simulation)
    f, axes = newfig(gridbox);
    plotsim(axes, sim)
    [Colorbar(f[i, 2], axes[i].scene.plots[2].plots[1]) for i in 1:length(axes)]
    f
end

f, axes = newfig(gridbox);
ax1, ax2, ax3 = axes
f = plotsim(sim)
plotsim(axes, sim)
plotsim(axes, sim.u, divergence, sim.p)

simcopy = copy(sim)
gp = griddedplot!(axes[1], sim.u)
pl = gp.plots[1]
gp = griddedplot!(ax1, uobs)
gp = griddedplot!(ax2, sim.v)
gp = griddedplot!(ax2, sim.p)
gp = griddedplot!(ax3, div)
contourf!(getxs(div), getys(div), div.mesh)

da = DimArray(u)
f = Figure()
ax = Axis(f[1,1])
bn.vertcut!(ax, u, [0, 1, 2, 4, 5, 8])
f

bn.vertcut(sim.u, [0, 1, 2, 4, 5, 8, 16], canal)
bn.vertcut(sim.u, [0, 2, 3.5, 4, 4.2, 8, 16], canal)
bn.vertcut(sim.v, [0, 1, 2, 4, 5, 8, 16], canal)
bn.vertcut(sim.v, [0, 2, 4, 10], canal)
bn.vertcut(sim.v, [0, 10], canal)
bn.vertcut(sim.u, [0, 16], canal)
bn.vertcut(VortGrid(sim.u, sim.v), [0, 1, 2, 4, 5, 8, 16], canal)
bn.vertcut(DivGrid(sim.u, sim.v), [0, 1, 2, 4, 4.2, 8, 16], canal)
bn.vertcut(sim.p, [0, 1, 2, 4, 5, 8, 16], canal)

bn.flowrate(sim.u[1, 2:end-1], sim.u.grid.dy)
plot(getxs(sim.u), bn.flowrate(sim.u))

plot(sim.u[end, :], getys(sim.u))
plot!(sim.u[end-1, :], getys(sim.u))
current_figure()
@testset "BottleNeckFlow.jl" begin
    # Write your tests here.
end