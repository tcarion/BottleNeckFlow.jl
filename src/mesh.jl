abstract type AbstractGrid{T} <: AbstractArray{T, 2} end

abstract type AbstractGridGhost{T} <: AbstractGrid{T} end
struct GridBox
    L::Float32
    D::Float32
    # Number of cells in x direction
    m::Int
    # Number of cells in y direction
    n::Int
    dx::Float64
    dy::Float64
end
function GridBox(L, D, m::Int, n::Int)
    dx = L / m
    dy = D / n
    GridBox(L, D, m, n, dx, dy)
end
GridBox(L, D, dx::AbstractFloat, dy::AbstractFloat) = GridBox(L, D, convert(Int, L / dx), convert(Int, D / dy))
GridBox(canal::CanalConfig, dx::AbstractFloat, dy::AbstractFloat) = GridBox(canal.L, canal.d, dx, dy)
GridBox(canal::CanalConfig, h::AbstractFloat) = GridBox(canal, h, h)
Base.size(gb::GridBox) = (gp.m, gb.n)
function Base.show(io::IO, gb::GridBox)
    println(io, "### GridBox:")
    for names in fieldnames(typeof(gb))
        println(io, "\t $(string(names)): $(getfield(gb, names))")
    end
end
from_ncf(::Type{GridBox}, fname::String) = _convert(GridBox, fname)

Base.parent(g::AbstractGrid) = parent(g.mesh)
Base.size(g::AbstractGrid) = size(parent(g))
Base.getindex(g::AbstractGrid, i::Int) = getindex(parent(g), i)
Base.getindex(g::AbstractGrid, i::Int, j::Int) = getindex(parent(g), i, j)
Base.setindex!(g::AbstractGrid, v,  i::Int) = setindex!(parent(g), v, i)
Base.setindex!(g::AbstractGrid, v,  i::Int, j::Int) = setindex!(parent(g), v, i, j)
# Base.copy(g::AbstractGridGhost) = typeof(g)(copy(parent(g)), g.grid, (copy(g.ghost[1]), (copy(g.ghost[2]))), g.ghostdim)
Base.copy(g::AbstractGrid) = typeof(g)(copy(parent(g)), g.grid, ghostinds1(g), ghostinds2(g))

getxs(g::AbstractGrid) = 0.5*g.grid.dx:g.grid.dx:g.grid.L-0.5*g.grid.dx
getys(g::AbstractGrid) = (-g.grid.D + g.grid.dx)*0.5:g.grid.dy:(g.grid.D - g.grid.dx)*0.5

ghostinds1(g::AbstractGridGhost) = g.ghostinds1
ghostinds2(g::AbstractGridGhost) = g.ghostinds2

average(v) = sum(v) / length(v)
# struct UGrid{T} <: AbstractGrid{T}
#     mesh::Matrix{T}
#     grid::GridBox
#     # ghost0::Vector{T}
#     # ghostn::Vector{T}
# end
struct UGrid{T} <: AbstractGridGhost{T}
    mesh::Matrix{T}
    grid::GridBox
    # ghost::Tuple{Vector{T}, Vector{T}}
    ghostinds1::AbstractVector{CartesianIndex}
    ghostinds2::AbstractVector{CartesianIndex}
end

# function UGrid{T}(gb::GridBox) where T <: AbstractFloat
#     UGrid(zeros(T, gb.m, gb.n), gb, zeros(T, gb.n), zeros(T, gb.n))
# end
# function UGrid{T}(gb::GridBox) where T <: AbstractFloat
#     UGrid(zeros(T, gb.m + 1, gb.n + 2), gb)
# end
function UGrid(matrix::Matrix{T}, gb::GridBox) where T <: AbstractFloat
    nx, ny = size(matrix)
    # ghost = (zeros(T, nx), zeros(T, nx))
    g1 = [CartesianIndex(i, 1) for i in 1:nx]
    g2 = [CartesianIndex(i, ny) for i in 1:nx]
    # UGrid{T}(matrix, gb, CartesianIndices((1:nx, 1)), CartesianIndices((1:nx, ny)))
    UGrid{T}(matrix, gb, g1, g2)
end
function UGrid{T}(gb::GridBox) where T <: AbstractFloat
    nx = gb.m + 1
    ny = gb.n + 2
    mesh = zeros(T, nx, ny)
    UGrid(mesh, gb)
end

UGrid(gb::GridBox) = UGrid{Float64}(gb)
atsw(::Type{<:UGrid}, i, j) = CartesianIndex(i, j-1)
atse(::Type{<:UGrid}, i, j) = CartesianIndex(i+1, j-1)
atne(::Type{<:UGrid}, i, j) = CartesianIndex(i+1, j)
atnw(::Type{<:UGrid}, i, j) = CartesianIndex(i, j)
atleft(::Type{<:UGrid}, i, j) = CartesianIndex(i-1, j-1)
atright(::Type{<:UGrid}, i, j) = CartesianIndex(i, j-1)
# projectionV(g::UGrid, i, j) = 0.25 * (g[i, j] + g[i, j-1] + g[i+1, j - 1] + g[i+1, j])
# projonV(::UGrid, i, j) = CartesianIndex.([(i, j-1), (i+1, j-1), (i+1, j), (i, j)])
getxs(g::UGrid) = 0:g.grid.dx:g.grid.L
getys(g::UGrid) = -g.grid.D*0.5 - 0.5*g.grid.dy: g.grid.dy :g.grid.D*0.5 + 0.5*g.grid.dy
# inrange(p::UGrid) = 2:size(p, 1)-1, 1:size(p, 2)
function trueinds(g::UGrid)
    nx, ny = size(g)
    CartesianIndices((nx, 2:ny-1))
end
function insideinds(g::UGrid) 
    nx, ny = size(g)
    CartesianIndices((2:nx-1, 2:ny-1))
end
# trueview(g::UGrid) = @view g.mesh[:, 2:end-1]

struct VGrid{T} <: AbstractGridGhost{T}
    mesh::Matrix{T}
    grid::GridBox
    # ghost::Tuple{Vector{T}, Vector{T}}
    # ghostdim::Int
    ghostinds1::AbstractVector{CartesianIndex}
    ghostinds2::AbstractVector{CartesianIndex}
end
function VGrid(matrix::Matrix{T}, gb::GridBox) where T <: AbstractFloat
    nx, ny = size(matrix)
    g1 = [CartesianIndex(1, j) for j in 1:ny]
    g2 = [CartesianIndex(nx, j) for j in 1:ny]
    # ghost = (zeros(T, ny), zeros(T, ny))
    VGrid{T}(matrix, gb, g1, g2)
end
function VGrid{T}(gb::GridBox) where T <: AbstractFloat
    nx = gb.m + 2
    ny = gb.n + 1
    mesh = zeros(T, nx, ny)
    VGrid(mesh, gb)
end
# struct VGrid{T} <: AbstractGrid{T}
#     mesh::Matrix{T}
#     grid::GridBox
# end
# function VGrid{T}(gb::GridBox) where T <: AbstractFloat
#     VGrid(zeros(T, gb.m + 2, gb.n + 1), gb)
# end
VGrid(gb::GridBox) = VGrid{Float64}(gb)
# projectionU(g::VGrid, i, j) = 0.25 * (g[i, j] + g[i, j+1] + g[i-1, j + 1] + g[i-1, j]) 
# projonU(::VGrid, i, j) = CartesianIndex.([(i-1, j), (i, j), (i, j+1), (i-1, j+1)])
atsw(::Type{<:VGrid}, i, j) = CartesianIndex(i-1, j)
atse(::Type{<:VGrid}, i, j) = CartesianIndex(i, j)
atne(::Type{<:VGrid}, i, j) = CartesianIndex(i, j+1)
atnw(::Type{<:VGrid}, i, j) = CartesianIndex(i-1, j+1)
attop(::Type{<:VGrid}, i, j) = CartesianIndex(i-1, j)
atbot(::Type{<:VGrid}, i, j) = CartesianIndex(i-1, j-1)

getxs(g::VGrid) = -g.grid.dx*0.5:g.grid.dx:g.grid.L+g.grid.dx*0.5
getys(g::VGrid) = -g.grid.D*0.5:g.grid.dy:g.grid.D*0.5

inrange(v::VGrid) = 1:size(v)[1], 2:size(v)[2]-1
function trueinds(g::VGrid)
    nx, ny = size(g)
    CartesianIndices((2:nx-1, ny))
end
function insideinds(g::VGrid) 
    nx, ny = size(g)
    CartesianIndices((2:nx-1, 2:ny-1))
end

neighbours(g::Type{<:AbstractGridGhost}, i, j) = [atsw(g, i, j), atse(g, i, j), atne(g, i, j), atnw(g, i, j)]
neighbours(g::Type{<:AbstractGridGhost},I::CartesianIndex) = neighbours(g, Tuple(I)...)
truegrid(g::AbstractGridGhost) = g[trueinds(g)]
atsw(g::Type{<:AbstractGridGhost}, I::CartesianIndex) = atsw(g, Tuple(I)...)
atse(g::Type{<:AbstractGridGhost}, I::CartesianIndex) = atse(g, Tuple(I)...)
atne(g::Type{<:AbstractGridGhost}, I::CartesianIndex) = atne(g, Tuple(I)...)
atnw(g::Type{<:AbstractGridGhost}, I::CartesianIndex) = atnw(g, Tuple(I)...)

# trueview(g::VGrid) = @view g.mesh[2:end-1, :]

struct PGrid{T} <: AbstractGrid{T}
    mesh::Matrix{T}
    grid::GridBox
end
# function PGrid(matrix::Matrix{T}, gb::GridBox) where T
#     PGrid{T}(matrix, gb)
# end
function PGrid{T}(gb::GridBox) where T <: AbstractFloat
    PGrid(zeros(T, gb.m, gb.n), gb)
end
PGrid(gb::GridBox) = PGrid{Float64}(gb)

# inrange(p::PGrid) = 2:size(p)[1]-1, 2:size(p)[2]-1
function insideinds(g::PGrid) 
    nx, ny = size(g)
    CartesianIndices((2:nx-1, 2:ny-1))
end
# getxs(g::PGrid) = 0+g.grid.dx*0.5:g.grid.dx:g.grid.L-g.grid.dx*0.5
# getys(g::PGrid) = -g.grid.D*0.5+0.5*g.grid.dy:g.grid.dy:g.grid.D*0.5-0.5*g.grid.dy

# trueview(g::PGrid) = @view g.mesh[:, :]

struct VortGrid{T} <: AbstractGrid{T}
    mesh::Matrix{T}
    grid::GridBox
end
function VortGrid{T}(gb::GridBox) where T <: AbstractFloat
    VortGrid(zeros(T, gb.m, gb.n), gb)
end
VortGrid(gb::GridBox) = VortGrid{Float64}(gb)
getxs(g::VortGrid) = 0+g.grid.dx*0.5:g.grid.dx:g.grid.L-g.grid.dx*0.5
getys(g::VortGrid) = -g.grid.D*0.5+0.5*g.grid.dy:g.grid.dy:g.grid.D*0.5-0.5*g.grid.dy

function VortGrid(u::UGrid, v::VGrid)
    vort = VortGrid(u.grid)
    for I in eachindex(vort)
        i,j = Tuple(I)
        vsleft = v[neighbours(UGrid, atleft(VortGrid, I))]
        vsright = v[neighbours(UGrid, atright(VortGrid, I))]
        usbot = u[neighbours(VGrid, atbot(VortGrid, I))]
        ustop = u[neighbours(VGrid, attop(VortGrid, I))]
        # dvdx = (projectionU(v, i+1, j) - projectionU(v, i, j)) / v.grid.dx
        dvdx = (average(vsright) - average(vsleft)) / v.grid.dx
        # dudy = (projectionV(u, i, j+1) - projectionV(u, i, j)) / v.grid.dy
        dudy = (average(ustop) - average(usbot)) / u.grid.dy
        vort[I] = dvdx - dudy
    end
    vort
end

struct DivGrid{T} <: AbstractGrid{T}
    mesh::Matrix{T}
    grid::GridBox
end
function DivGrid{T}(gb::GridBox) where T <: AbstractFloat
    DivGrid(zeros(T, gb.m, gb.n), gb)
end
DivGrid(gb::GridBox) = DivGrid{Float64}(gb)
getxs(g::DivGrid) = 0+g.grid.dx*0.5:g.grid.dx:g.grid.L-g.grid.dx*0.5
getys(g::DivGrid) = -g.grid.D*0.5+0.5*g.grid.dy:g.grid.dy:g.grid.D*0.5-0.5*g.grid.dy

function DivGrid(u::UGrid, v::VGrid)
    div = DivGrid(u.grid)
    for ei in eachindex(div)
        i,j = Tuple(ei)
        dudx = (u[i+1, j] - u[i, j]) / v.grid.dx
        dvdy = (v[i, j+1] - v[i, j]) / v.grid.dy
        div[ei] = dudx + dvdy
    end
    div
end

atleft(::Type{<:AbstractGrid}, i, j) = CartesianIndex(i, j+1)
atbot(::Type{<:AbstractGrid}, i, j) = CartesianIndex(i+1, j)
atright(::Type{<:AbstractGrid}, i, j) = CartesianIndex(i+1, j+1)
attop(::Type{<:AbstractGrid}, i, j) = CartesianIndex(i+1, j+1)

atleft(g, I::CartesianIndex) = atleft(g, Tuple(I)...)
atbot(g, I::CartesianIndex) = atbot(g, Tuple(I)...)
atright(g, I::CartesianIndex) = atright(g, Tuple(I)...)
attop(g, I::CartesianIndex) = attop(g, Tuple(I)...)
# trueview(g::VortGrid) = @view g.mesh[:, :]

function DimensionalData.DimArray(g::AbstractGrid)
    DimArray(parent(g), (X(getxs(g) |> collect), Y(getys(g) |> collect)))
end