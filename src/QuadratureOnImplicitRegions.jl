module QuadratureOnImplicitRegions


using FastGaussQuadrature,LinearAlgebra,ForwardDiff

import FiniteDiff.finite_difference_gradient! #non-allocating gradient

export algoim_nodes_weights
# export algoim_quad #nneds work

include("algoim.jl")

# include("integrate_f.jl") #Needs some work
include("algoim_functions.jl")


end
