# QuadratureOnImplicitRegions

[![Build Status](https://github.com/hmegh/QuadratureOnImplicitRegions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hmegh/QuadratureOnImplicitRegions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://hmegh.github.io/QuadratureOnImplicitRegions.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://hmegh.github.io/QuadratureOnImplicitRegions.jl/dev/)

This package implements a quadrature method on implicitly defined regions following the algorithm in: 

[R. I. Saye, High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles, SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015).](https://epubs.siam.org/doi/10.1137/140966290).


---

## Examples: 

Let $\Omega=[0,1]^2$ and let $\psi(x,y)=x^2+y^2-1$. Our goal is to create quadrature nodes and weights on the subdomains: 

$$\Omega_1=\left\\{(x,y)\in \Omega : \psi(x,y)<0\right\\},\qquad 
\Omega_2=\left\\{(x,y)\in \Omega : \psi(x,y)>0\right\\}.$$

```julia
using QuadratureOnImplicitRegions

ψ(x)= x'*x-1.0 
a,b=zeros(2), ones(2) #the unit interval. 
quad_order=10

#the nodes and weights on Ω₁
xy1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)
#the nodes and weights on Ω₂
xy2,w2=algoim_nodes_weights(ψ,+1.0, a,b,quad_order)
```
To plot the nodes, please see [this tutorial](https://github.com/Hmegh/QuadratureOnImplicitRegions.jl/blob/main/tutorial/circle_and_sphere.jl).

<p align="center">
  <img src="https://github.com/Hmegh/QuadratureOnImplicitRegions.jl/assets/8241188/8926d082-3b1c-48cb-a888-3882b1288f7f" width="250" 
     height=auto/>
</p>

The same syntax can be used for higher dimensional regions. For example, in the case of the intersection of the unit sphere and unit cube, we only need to adjust `a` and `b`:

```julia
using QuadratureOnImplicitRegions

ψ(x)= x'*x-1.0 
a,b=zeros(3), ones(3) #the unit cube. 
quad_order=5 

xyz1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)
```
For the outer region, we only need to change `-1.0` to `1.0`


<p align="center">
  <img src="https://github.com/Hmegh/QuadratureOnImplicitRegions.jl/assets/8241188/43354dab-7818-46eb-8ee2-c65b394b0369" width="250" 
     height=auto/>
</p>


---

## Roadmap 

In the future, I plan on adding the following features/improvements. Suggestions and/or contributions are appreciated. 

- [ ] Add support for regions defined as $\psi(\mathbf{x})=0$ (e.g. sphere in $\mathbb{R}^3$). 
- [ ] Performance improvement: See the [v0.2.0 release notes](https://github.com/HMegh/QuadratureOnImplicitRegions.jl/releases/tag/v0.2.0) for a comparison with the previous release. This is still far from optimal. 
- [ ] Add the option to integrate a function $f$ without the need to store the integration nodes. This is useful when the nodes are not re-used. 
- [ ] Adaptivity? It would be cool to have a [QuadGK](https://juliamath.github.io/QuadGK.jl/stable/)-like function, where the error is estimated using a Gauss-Konrod rule. 
