module BottleNeckFlow

# using RecipesBase
using Makie
using Logging
using NCDatasets
using Dates
# using DataStructures
# using Rasters

export
    Simulation, 
    CanalConfig, δ, SimParam,
    GridBox,
    AbstractGrid,
    UGrid,
    VGrid,
    PGrid,
    VortGrid,
    DivGrid,
    getxs,
    getys,
    poiseuille,
    poiseuille!,
    fillghost!,
    noslip,
    naturalbound,
    inbump,
    step_euler!,
    gauss_seidel!,
    residual,
    step_poisson!,
    poisson_project!,
    runsim!,
    create,
    add_to_nc

include("utils.jl")
include("constants.jl")
include("canal.jl")
include("mesh.jl")
include("simulation.jl")
include("poiseuille.jl")
include("schemes.jl")
include("diagnostic.jl")
include("io.jl")
include("run.jl")

include("recipes.jl")

end
