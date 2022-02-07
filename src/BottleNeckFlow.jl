module BottleNeckFlow

# using RecipesBase
using Makie

export CanalConfig, δ,
    GridBox,
    AbstractGrid,
    UGrid,
    VGrid,
    getxs,
    getys,
    poiseuille,
    poiseuille!

include("constants.jl")
include("canal.jl")
include("mesh.jl")
include("poiseuille.jl")
include("schemes.jl")

include("recipes.jl")

end
