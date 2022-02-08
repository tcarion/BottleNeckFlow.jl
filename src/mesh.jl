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

abstract type AbstractGrid{T} <: AbstractArray{T, 2} end

abstract type AbstractGridGhost{T} <: AbstractGrid{T} end
Base.size(g::AbstractGrid) = size(g.mesh)
Base.getindex(g::AbstractGrid, i::Int) = getindex(g.mesh, i)
Base.getindex(g::AbstractGrid, i::Int, j::Int) = getindex(g.mesh, i, j)
Base.setindex!(g::AbstractGrid, v,  i::Int) = setindex!(g.mesh, v, i)
Base.setindex!(g::AbstractGrid, v,  i::Int, j::Int) = setindex!(g.mesh, v, i, j)

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
function UGrid{T}(gb::GridBox) where T <: AbstractFloat
    nx = gb.m + 1
    ny = gb.n
    mesh = zeros(T, nx, ny)
    ghost = (zeros(T, nx), zeros(T, nx))
    UGrid(mesh, gb, ghost, 2)
end

UGrid(gb::GridBox) = UGrid{Float64}(gb)
average2d(g::UGrid, i, j) = 0.25 * (g[i, j] + g[i, j-1] + g[i+1, j - 1] + g[i+1, j]) 
# getxs(g::UGrid) = 0:g.grid.dx:g.grid.L
# getys(g::UGrid) = -g.grid.D*0.5 + 0.5*g.grid.dy:g.grid.dy:g.grid.D*0.5 - 0.5*g.grid.dy

# trueview(g::UGrid) = @view g.mesh[:, 2:end-1]

struct VGrid{T} <: AbstractGridGhost{T}
    mesh::Matrix{T}
    grid::GridBox
    ghost::Tuple{Vector{T}, Vector{T}}
    ghostdim::Int
end
function VGrid{T}(gb::GridBox) where T <: AbstractFloat
    nx = gb.m
    ny = gb.n + 1
    mesh = zeros(T, nx, ny)
    ghost = (zeros(T, ny), zeros(T, ny))
    VGrid(mesh, gb, ghost, 1)
end
# struct VGrid{T} <: AbstractGrid{T}
#     mesh::Matrix{T}
#     grid::GridBox
# end
# function VGrid{T}(gb::GridBox) where T <: AbstractFloat
#     VGrid(zeros(T, gb.m + 2, gb.n + 1), gb)
# end
VGrid(gb::GridBox) = VGrid{Float64}(gb)
average2d(g::VGrid, i, j) = 0.25 * (g[i, j] + g[i, j+1] + g[i-1, j + 1] + g[i-1, j]) 

# getxs(g::VGrid) = 0+g.grid.dx*0.5:g.grid.dx:g.grid.L-g.grid.dx*0.5
# getys(g::VGrid) = -g.grid.D*0.5:g.grid.dy:g.grid.D*0.5

# trueview(g::VGrid) = @view g.mesh[2:end-1, :]

struct PGrid{T} <: AbstractGrid{T}
    mesh::Matrix{T}
    grid::GridBox
end
function PGrid{T}(gb::GridBox) where T <: AbstractFloat
    PGrid(zeros(T, gb.m, gb.n), gb)
end
PGrid(gb::GridBox) = PGrid{Float64}(gb)
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
        dvdx = (average2d(v, i+1, j) - average2d(v, i, j)) / v.grid.dx
        dudy = (average2d(u, i, j+1) - average2d(u, i, j)) / v.grid.dy
        vort[ei] = dvdx - dudy
    end
    vort
end

# trueview(g::VortGrid) = @view g.mesh[:, :]