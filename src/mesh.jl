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

abstract type AbstractGrid{T} <: AbstractArray{T, 2} end

abstract type AbstractGridGhost{T} <: AbstractGrid{T} end
Base.size(g::AbstractGrid) = size(g.mesh)
Base.getindex(g::AbstractGrid, i::Int) = getindex(g.mesh, i)
Base.getindex(g::AbstractGrid, i::Int, j::Int) = getindex(g.mesh, i, j)
Base.setindex!(g::AbstractGrid, v,  i::Int) = setindex!(g.mesh, v, i)
Base.setindex!(g::AbstractGrid, v,  i::Int, j::Int) = setindex!(g.mesh, v, i, j)
Base.copy(g::AbstractGridGhost) = typeof(g)(copy(g.mesh), g.grid, (copy(g.ghost[1]), (copy(g.ghost[2]))), g.ghostdim)
Base.copy(g::AbstractGrid) = typeof(g)(copy(g.mesh), g.grid)

getxs(g::AbstractGrid) = 0.5*g.grid.dx:g.grid.dx:g.grid.L-0.5*g.grid.dx
getys(g::AbstractGrid) = (-g.grid.D + g.grid.dx)*0.5:g.grid.dy:(g.grid.D - g.grid.dx)*0.5

function Base.getindex(g::AbstractGridGhost, i::Int, j::Int)
    try
        getindex(g.mesh, i, j)
    catch e
        if e isa BoundsError
            gd = g.ghostdim
            od = gd == 1 ? 2 : 1
            if e.i[g.ghostdim] == 0
                getindex(g.ghost[1], e.i[od])
            elseif e.i[g.ghostdim] == size(g)[gd] + 1
                getindex(g.ghost[2], e.i[od])
            else
                rethrow(e)
            end
        else
            rethrow(e)
        end
    end
end
function Base.setindex!(g::AbstractGridGhost, v, i::Int, j::Int)
    try
        setindex!(g.mesh, v, i, j)
    catch e
        if e isa BoundsError
            gd = g.ghostdim
            od = gd == 1 ? 2 : 1
            if e.i[g.ghostdim] == 0
                setindex!(g.ghost[1], v, e.i[od])
            elseif e.i[g.ghostdim] == size(g)[gd] + 1
                setindex!(g.ghost[2], v, e.i[od])
            else
                rethrow(e)
            end
        else
            rethrow(e)
        end
    end
end
# struct UGrid{T} <: AbstractGrid{T}
#     mesh::Matrix{T}
#     grid::GridBox
#     # ghost0::Vector{T}
#     # ghostn::Vector{T}
# end
struct UGrid{T} <: AbstractGridGhost{T}
    mesh::Matrix{T}
    grid::GridBox
    ghost::Tuple{Vector{T}, Vector{T}}
    ghostdim::Int
end

# function UGrid{T}(gb::GridBox) where T <: AbstractFloat
#     UGrid(zeros(T, gb.m, gb.n), gb, zeros(T, gb.n), zeros(T, gb.n))
# end
# function UGrid{T}(gb::GridBox) where T <: AbstractFloat
#     UGrid(zeros(T, gb.m + 1, gb.n + 2), gb)
# end
function UGrid(matrix::Matrix{T}, gb::GridBox) where T <: AbstractFloat
    nx = size(matrix, 1)
    ghost = (zeros(T, nx), zeros(T, nx))
    UGrid{T}(matrix, gb, ghost, 2)
end
function UGrid{T}(gb::GridBox) where T <: AbstractFloat
    nx = gb.m + 1
    ny = gb.n
    mesh = zeros(T, nx, ny)
    UGrid(mesh, gb)
end

UGrid(gb::GridBox) = UGrid{Float64}(gb)
projectionV(g::UGrid, i, j) = 0.25 * (g[i, j] + g[i, j-1] + g[i+1, j - 1] + g[i+1, j])
getxs(g::UGrid) = 0:g.grid.dx:g.grid.L
getys(g::UGrid) = -g.grid.D*0.5 + 0.5*g.grid.dy:g.grid.dy:g.grid.D*0.5 - 0.5*g.grid.dy

inrange(p::UGrid) = 2:size(p)[1]-1, 1:size(p)[2]

# trueview(g::UGrid) = @view g.mesh[:, 2:end-1]

struct VGrid{T} <: AbstractGridGhost{T}
    mesh::Matrix{T}
    grid::GridBox
    ghost::Tuple{Vector{T}, Vector{T}}
    ghostdim::Int
end
function VGrid(matrix::Matrix{T}, gb::GridBox) where T <: AbstractFloat
    ny = size(matrix, 2)
    ghost = (zeros(T, ny), zeros(T, ny))
    VGrid{T}(matrix, gb, ghost, 1)
end
function VGrid{T}(gb::GridBox) where T <: AbstractFloat
    nx = gb.m
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
projectionU(g::VGrid, i, j) = 0.25 * (g[i, j] + g[i, j+1] + g[i-1, j + 1] + g[i-1, j]) 

getxs(g::VGrid) = 0+g.grid.dx*0.5:g.grid.dx:g.grid.L-g.grid.dx*0.5
getys(g::VGrid) = -g.grid.D*0.5:g.grid.dy:g.grid.D*0.5

inrange(v::VGrid) = 1:size(v)[1], 2:size(v)[2]-1


# trueview(g::VGrid) = @view g.mesh[2:end-1, :]

struct PGrid{T} <: AbstractGrid{T}
    mesh::Matrix{T}
    grid::GridBox
end
function PGrid(matrix::Matrix{T}, gb::GridBox) where T
    PGrid{T}(matrix, gb)
end
function PGrid{T}(gb::GridBox) where T <: AbstractFloat
    PGrid(zeros(T, gb.m, gb.n), gb)
end
PGrid(gb::GridBox) = PGrid{Float64}(gb)

inrange(p::PGrid) = 2:size(p)[1]-1, 2:size(p)[2]-1
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
    for ei in eachindex(vort)
        i,j = Tuple(ei)
        dvdx = (projectionU(v, i+1, j) - projectionU(v, i, j)) / v.grid.dx
        dudy = (projectionV(u, i, j+1) - projectionV(u, i, j)) / v.grid.dy
        vort[ei] = dvdx - dudy
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

# trueview(g::VortGrid) = @view g.mesh[:, :]