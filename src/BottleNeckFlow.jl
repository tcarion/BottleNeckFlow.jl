module BottleNeckFlow

# using RecipesBase
using Makie
using Logging

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
    runsim!

include("constants.jl")
include("canal.jl")
include("mesh.jl")
include("simulation.jl")
include("poiseuille.jl")
include("schemes.jl")
include("diagnostic.jl")
include("run.jl")

include("recipes.jl")

end
