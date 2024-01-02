using QuadratureOnImplicitRegions
using Test

@testset "QuadratureOnImplicitRegions.jl" begin
    # Write your tests here.
    ψ(x)=x'*x-1

    #quarter of a circle
    a,b=zeros(2), ones(2)
    x,w=algoim_nodes_weights(ψ,-1.0, a,b,10)
    @test sum(w) ≈ π/4

    #1/8 of a sphere
    a,b=zeros(3), ones(3)
    x,w=algoim_nodes_weights(ψ,-1.0, a,b,10)
    @test sum(w) ≈ π/6
end
