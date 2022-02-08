using BottleNeckFlow
using CairoMakie
using GLMakie
using Test

bn = BottleNeckFlow

H = 1.
canal = CanalConfig(H)
params = SimParam()
Re = 250
ubar = 3.
nu = ubar * canal.L / Re
rho = 1000
mu = nu * rho

canobs = Observable(CanalConfig(2.))
x = range(0, canal.L, length = 101)
xobs = range(0, L, length = 50) |> Observable
y = δ(x, canal)

gridbox = GridBox(canal, .5)
u = UGrid(gridbox)
xu = getxs(u) |> collect
yu = getys(u) |> collect

v = VGrid(gridbox)
getxs(v) |> collect
getys(v) |> collect

p = PGrid(gridbox)
getxs(p) |> collect
getys(p) |> collect

poiseuille!(u, canal, params)

vort = VortGrid(u.grid)
vort = VortGrid(u, v)

sim = Simulation(params, gridbox)
poiseuille!(sim)

canalplot(x, canal)
f = Figure()
ax = f[1, 1] = Axis(f, aspect = DataAspect())
canalplot!(ax, x, canal)
canalplot!(ax, xobs, canobs)

canobs[] = CanalConfig(1.5)
xobs[] = range(0, L-2, length= 20)

contourf!(ax, xu, yu, u, colormap = :jet)
contourf(xu, yu, u, colormap = :jet)
cont = contourf!(ax, getxs(vort), getys(vort), vort .|> abs, colormap = :jet)
cont = contourf(getxs(vort), getys(vort), vort .|> abs, colormap = :jet)
Colorbar(f[1, 2], cont)
xlims!(ax, (0, u.grid.L))
ylims!(ax, (0, u.grid.D))

gp = griddedplot!(ax, u)
@testset "BottleNeckFlow.jl" begin
    # Write your tests here.
end