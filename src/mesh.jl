struct GridBox
    L::Float32
    H::Float32
    m::Int
    n::Int
    dx::Float64
    dy::Float64
end
function GridBox(L, H, m::Int, n::Int)
    dx = L / (m - 1)
    dy = H / (n - 1)
    GridBox(L, H, m, n, dx, dy)
end
GridBox(L, H, dx::AbstractFloat, dy::AbstractFloat) = GridBox(L, H, L / dx + 1, H / dy + 1, dx, dy)

abstract type AbstractGrid{T} <: AbstractArray{T, 2} end

Base.size(g::AbstractGrid) = size(g.mesh)
Base.getindex(g::AbstractGrid, i::Int) = getindex(g.mesh, i)
Base.getindex(g::AbstractGrid, i::Int, j::Int) = getindex(g.mesh, i, j)
Base.setindex!(g::AbstractGrid, v,  i::Int) = setindex!(g.mesh, v, i)
Base.setindex!(g::AbstractGrid, v,  i::Int, j::Int) = setindex!(g.mesh, v, i, j)


struct UGrid{T} <: AbstractGrid{T}
    mesh::Matrix{T}
    grid::GridBox
    # ghost0::Vector{T}
    # ghostn::Vector{T}
end
# function UGrid{T}(gb::GridBox) where T <: AbstractFloat
#     UGrid(zeros(T, gb.m, gb.n), gb, zeros(T, gb.n), zeros(T, gb.n))
# end
function UGrid{T}(gb::GridBox) where T <: AbstractFloat
    UGrid(zeros(T, gb.m + 1, gb.n), gb)
end
UGrid(gb::GridBox) = UGrid{Float64}(gb)
getxs(g::UGrid) = g.grid.dx/2:g.grid.dx:g.grid.L 
getys(g::UGrid) = 0:g.grid.dy:g.grid.H

innerview(g::UGrid) = @view g.mesh[2:end-1, :]


struct VGrid{T} <: AbstractGrid{T}
    mesh::Matrix{T}
    grid::GridBox
end
function VGrid{T}(gb::GridBox) where T <: AbstractFloat
    VGrid(zeros(T, gb.m, gb.n + 1), gb)
end
VGrid(gb::GridBox) = VGrid{Float64}(gb)
getxs(g::VGrid) = 0:g.grid.dx:g.grid.L 
getys(g::VGrid) = g.grid.dy/2:g.grid.dy:g.grid.H

innerview(g::VGrid) = @view g.mesh[:, 2:end-1]