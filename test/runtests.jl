using QuadratureOnImplicitRegions
using Test
using SpecialFunctions


@testset "QuadratureOnImplicitRegions.jl" begin
    # Write your tests here.
    ψ(x)=x'*x-1


    ## Testing the weights (should sum up to the area): 

    #quarter of a circle
    a,b=zeros(2), ones(2)
    x,w=algoim_nodes_weights(ψ,-1.0, a,b,10)
    @test sum(w) ≈ π/4


    #ellipse of radii 1,1/2
    ell(x) = x[1]^2+(2*x[2])^2 -1.0
    a,b=-2ones(2),2ones(2)
    x,w=algoim_nodes_weights(ell,-1.0, a,b,10)
    @test sum(w) ≈ π/2

    #1/8 of a sphere
    a,b=zeros(3), ones(3)
    x,w=algoim_nodes_weights(ψ,-1.0, a,b,10)
    @test sum(w) ≈ π/6


    ## Testing the quadrature rule for polynomials
    #This help to make sure that the nodes are mapped approperiately. 
    #Note: Keep in mind that this quadrature is not exact (even for polynomials). 

    #quarter of a circle
    a,b=zeros(2), ones(2)
    x,w=algoim_nodes_weights(ψ,-1.0, a,b,11)

    # ∫∫ x^i dxdy on the quarter circle 
    #   = √π/4*Γ((i+1)/2)/Γ(2+i/i)
    for i=1:10
        exact_integral= √π/4 * gamma((i+1)/2)/gamma(2+i/2)
        approximate_integral=  w'x[1,:].^i 
        @test exact_integral≈ approximate_integral
    end

    # ∫∫ x^i y^j dxdy on the quarter circle 
    #   = Γ((i+1)/2)Γ((j+3)/2)/(2(j+1)Γ(2+(i+j)/2))

    for i=1:10
        for j=1:10
            exact_integral=gamma((i+1)/2)*gamma((j+3)/2)/(2(j+1)*gamma(2+(i+j)/2))
            approximate_integral=w'*(x[1,:].^i .* x[2,:].^j)
            @test exact_integral≈ approximate_integral
        end
    end

end
