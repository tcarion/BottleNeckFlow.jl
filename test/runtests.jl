using BottleNeckFlow
using CairoMakie
using GLMakie
using Test

bn = BottleNeckFlow

H = 1.
L = 10*H
Re = 2000
ubar = 3.
nu = ubar * L / Re
rho = 1000
mu = nu * rho

canal = CanalConfig(H)
canobs = Observable(CanalConfig(2.))
x = range(0, L, length = 101)
xobs = range(0, L, length = 50) |> Observable
y = δ(x, canal)

gridbox = GridBox(L, H, 21, 11)
u = UGrid(gridbox)
xu = getxs(u) |> collect
yu = getys(u) |> collect

v = VGrid(gridbox)
getxs(v) |> collect
getys(v) |> collect


poiseuille!(u, canal, mu, -1.)


f = Figure()
ax = f[1, 1] = Axis(f, aspect = DataAspect())
canalplot!(ax, x, canal)
canalplot!(ax, xobs, canobs)

canobs[] = CanalConfig(1.5)
xobs[] = range(0, L-2, length= 20)

contourf!(ax, xu, yu, innerview(u), colormap = :jet)
xlims!(ax, (0, u.grid.L))
ylims!(ax, (0, u.grid.H))

gp = griddedplot!(ax, u)
@testset "BottleNeckFlow.jl" begin
    # Write your tests here.
end