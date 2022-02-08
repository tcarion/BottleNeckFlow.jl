module BottleNeckFlow

# using RecipesBase
using Makie

export
    Simulation, 
    CanalConfig, δ, SimParam,
    GridBox,
    AbstractGrid,
    UGrid,
    VGrid,
    PGrid,
    VortGrid,
    getxs,
    getys,
    poiseuille,
    poiseuille!

include("constants.jl")
include("canal.jl")
include("mesh.jl")
include("simulation.jl")
include("poiseuille.jl")
include("schemes.jl")

include("recipes.jl")

end
