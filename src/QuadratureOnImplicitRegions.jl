module QuadratureOnImplicitRegions


using FastGaussQuadrature


export algoim_nodes_weights
export algoim_quad

include("algoim.jl")

include("integrate_f.jl")
include("algoim_functions.jl")


end
