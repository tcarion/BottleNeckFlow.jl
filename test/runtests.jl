using BottleNeckFlow
using CairoMakie
using Test

H = 1.
L = 10*H
canal = CanalConfig(H)

x = range(0, L, length = 101)

y = δ(x, canal)

f = Figure()
ax = f[1, 1] = Axis(f, aspect = DataAspect())
canalplot!(ax, x, canal)
@testset "BottleNeckFlow.jl" begin
    # Write your tests here.
end

lines!(ax, x, y)

lines!(ax, x, canal.H .- y)