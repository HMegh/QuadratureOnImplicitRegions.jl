var documenterSearchIndex = {"docs":
[{"location":"examples/","page":"Examples","title":"Examples","text":"The following are multiple examples showcasing the use of the main function algoim_nodes_weights. ","category":"page"},{"location":"examples/#2D-example","page":"Examples","title":"2D example","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this example, we consider the regions","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Omega_1=(xy)in 01^2mid  x^2+y^21qquad \nOmega_2= (xy)in 01^2mid  x^2+y^21","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Both regions are described by the level-set function psi(mathbfx)=Vert mathbfxVert^2-1. More precisely, ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Omega_1=(xy)in 01^2 mid psi(mathbfx)0qquad Omega_1=(xy)in 01^2 mid psi(mathbfx)0","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now, to set up the problem, we define psi, the domain 01times 01=a_1b_1times a_2b_2 and choose an order for the quadrature rule. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"ψ(x)= x'*x-1.0 \na,b=zeros(2), ones(2) #The domain [0,1] x [0,1]\nquad_order=10","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Next, we use algoim_nodes_weights to compute the quadrature rule on each region separately. Notice that the second argument represents the sign of psi on each region. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using QuadratureOnImplicitRegions\n\nxy1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)\nxy2,w2=algoim_nodes_weights(ψ,+1.0, a,b,quad_order)\n","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The first output xy1 is a matrix, where each column represents a quadrature node. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"xy1","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The second output w1 is a vector containing the weights. For instance, in this example, we expect the sum of w1 to be close to fracpi4. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"println(\"The sum of weights is: $(sum(w1)).\nThe error in the approximation of the area is $(pi/4- sum(w1)).\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To compute an integral on a region like ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"int_Omega_1 f(xy) dxdy","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Here’s a revised version of your sentence:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Evaluate the function f at the quadrature nodes, multiply by the corresponding weights, and then sum the results. This follows the same procedure as with the traditional Gaussian quadrature","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"f(x,y)=x*y #An example of a an integrand \napprox_integral= sum(f(xy1[1,i],xy1[2,i])*w1[i] for i in eachindex(w1))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"You can verify that the exact integral is frac18. Thus, the error is ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"approx_integral-1/8","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"You can plot the nodes using your favorite plotting library. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using CairoMakie\n\nfig=Figure()\nax = Axis(fig[1, 1], xlabel = \"x\", ylabel = \"y\",\n    title = \"Quadrature nodes\")\n\nsc1=scatter!(ax,xy1[1,:],xy1[2,:],color=:blue)\nsc2=scatter!(ax, xy2[1,:],xy2[2,:],color=:red)\n\nLegend(fig[1,2],[sc1, sc2],[\"nodes on Ω₁\",\"nodes on Ω₂\"])\n\nfig","category":"page"},{"location":"examples/#3D-example","page":"Examples","title":"3D example","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The same syntax can be used for regions in three dimensions (and even higher), as the following example shows. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using QuadratureOnImplicitRegions, CairoMakie\n\nψ(x)= x'*x-1.0 \na,b=zeros(3), ones(3) #the unit cube. \nquad_order=5 \n\nxyz1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)\n\nfig=Figure()\nax = Axis3(fig[1, 1], xlabel = \"x\", ylabel = \"y\",zlabel=\"z\",\n    title = \"Quadrature nodes on Ω₁\",\n    azimuth = -.2*pi,\n    elevation= pi/6 )\n\nscatter!(ax,xyz1[1,:],xyz1[2,:],xyz1[3,:],color=:blue)\n\nfig","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules =[QuadratureOnImplicitRegions]","category":"page"},{"location":"reference/#QuadratureOnImplicitRegions.d1_case-NTuple{5, Any}","page":"Reference","title":"QuadratureOnImplicitRegions.d1_case","text":"d1_case(ψ_list,s_list,domain,x_ref,w_ref)\n\nReturns a quadrature rule on the subset where si*ψi<0 (see the paper).\n\n\n\n\n\n","category":"method"},{"location":"reference/#QuadratureOnImplicitRegions.d1_count_subintervals-Tuple{Any, Any, Any}","page":"Reference","title":"QuadratureOnImplicitRegions.d1_count_subintervals","text":"d1_count_subintervals(ψ_list,s_list,domain)\n\nCounts how many sub-intervals of domain satisfiy siψi≥0\n\n\n\n\n\n","category":"method"},{"location":"reference/#QuadratureOnImplicitRegions.find_root-Union{Tuple{T}, Tuple{F}, Tuple{F, T, T}} where {F, T<:AbstractFloat}","page":"Reference","title":"QuadratureOnImplicitRegions.find_root","text":"find_root(ψ::F,a::T,b::T)\n\nFinds a root of a function ψ in the open interval (a,b)\n\n\n\n\n\n","category":"method"},{"location":"reference/#QuadratureOnImplicitRegions.find_roots-Union{Tuple{T}, Tuple{F}, Tuple{F, T, T}} where {F, T}","page":"Reference","title":"QuadratureOnImplicitRegions.find_roots","text":"find_roots(ψ::F,a::T,b::T) where {F,T}\n\nReturns multiple roots of ψ in the open interval [a,b]\n\n\n\n\n\n","category":"method"},{"location":"reference/#QuadratureOnImplicitRegions.find_roots-Union{Tuple{T}, Tuple{F}, Tuple{Vector{F}, T, T}} where {F, T}","page":"Reference","title":"QuadratureOnImplicitRegions.find_roots","text":"find_roots(ψ_list::Vector{F},a::T,b::T) where {F,T}\n\nReturns multiple roots of multiple functions ψ in the closed interval [a,b], it includes the endpoints by default even if they are not roots. \n\n\n\n\n\n","category":"method"},{"location":"reference/#QuadratureOnImplicitRegions.shifted_gl!-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, T, T, Vector{T}, Vector{T}}} where T","page":"Reference","title":"QuadratureOnImplicitRegions.shifted_gl!","text":"shifted_gl!(x_ref::Vector{T},w_ref::Vector{T},a::T,b::T,x::Vector{T},w::Vector{T})\n\nShifts the Gauss-Legendre rule (xref,wref) from (-1,1) to the interval [a,b] and stores them in x,w.\n\n\n\n\n\n","category":"method"},{"location":"reference/#QuadratureOnImplicitRegions.shifted_gl-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, T, T}} where T","page":"Reference","title":"QuadratureOnImplicitRegions.shifted_gl","text":"shifted_gl(x_ref::Vector{T},w_ref::Vector{T},a::T,b::T)\n\nShifts the Gauss-Legendre rule (xref,wref) from (-1,1) to the interval [a,b]\n\n\n\n\n\n","category":"method"},{"location":"reference/#QuadratureOnImplicitRegions.tensor_GL_rule-Union{Tuple{T}, Tuple{Matrix{T}, Vector{T}, Vector{T}}} where T","page":"Reference","title":"QuadratureOnImplicitRegions.tensor_GL_rule","text":"tensor_GL_rule(U::Matrix{T},x_ref::Vector{T},w_ref::Vector{T})\n\nReturns a matrix x and a vector w (multidimensional quadrature)     \n\n\n\n\n\n","category":"method"},{"location":"#QuadratureOnImplicitRegions.jl","page":"Home","title":"QuadratureOnImplicitRegions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuadratureOnImplicitRegions.jl is a package for generating high order quadrature rules on implicitly defined regions in hyperrectangles. ","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Simply add this package using the repl","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add QuadratureOnImplicitRegions","category":"page"},{"location":"#Tutorial-examples","page":"Home","title":"Tutorial examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See the examples section for a few examples","category":"page"},{"location":"#Reproducibility","page":"Home","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<details><summary>This documentation was built using these direct dependencies ,</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"}]
}
