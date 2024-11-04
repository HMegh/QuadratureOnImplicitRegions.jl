


The following are multiple examples showcasing the use of the main function `algoim_nodes_weights`. 



## 2D example

In this example, we consider the regions

$$\Omega_1=\{(x,y)\in [0,1]^2\mid  x^2+y^2<1\},\qquad 
\Omega_2=\{ (x,y)\in [0,1]^2\mid  x^2+y^2>1\}.$$

Both regions are described by the __level-set__ function $\psi(\mathbf{x})=\Vert \mathbf{x}\Vert^2-1$. More precisely, 

$$\Omega_1=\{(x,y)\in [0,1]^2\psi(\mathbf{x})<0\},\qquad \Omega_1=\{(x,y)\in [0,1]^2\psi(\mathbf{x})<0\}$$

Now, to set up the problem, we define $\psi$, the domain $[0,1]\times [0,1]=[a_1,b_1]\times [a_2,b_2]$ and choose an order for the quadrature rule. 

```@example simple
ψ(x)= x'*x-1.0 
a,b=zeros(2), ones(2) #The domain [0,1] x [0,1]
quad_order=10
```

Next, we use `algoim_nodes_weights` to compute the quadrature rule on each region separately. Notice that the second argument represents the sign of $\psi$ on each region. 

```@example simple
using QuadratureOnImplicitRegions

xy1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)
xy2,w2=algoim_nodes_weights(ψ,+1.0, a,b,quad_order)

```
The first output `xy1` is a matrix, where each column represents a quadrature node. 


```@example simple
xy1
```

The second output `w1` is a vector containing the weights. For instance, in this example, we expect the sum of `w1` to be close to $\frac{\pi}{4}$. 

```@example simple
println("The sum of weights is: $(sum(w1)).
The error in the approximation of the area is $(pi/4- sum(w1)).")
```
To compute an integral on a region like 

$$\int_{\Omega_1} f(x,y) dxdy,$$ 

Here’s a revised version of your sentence:

Evaluate the function $f$ at the quadrature nodes, multiply by the corresponding weights, and then sum the results. This follows the same procedure as with the traditional [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature)

```@example simple
f(x,y)=x*y #An example of a an integrand 
approx_integral= sum(f(xy1[1,i],xy1[2,i])*w1[i] for i in eachindex(w1))
```
You can verify that the exact integral is $\frac{1}{8}$. Thus, the error is 

```@example simple
approx_integral-1/8
```


You can plot the nodes using your favorite plotting library. 

```@example simple
using CairoMakie

fig=Figure()
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y",
    title = "Quadrature nodes")

sc1=scatter!(ax,xy1[1,:],xy1[2,:],color=:blue)
sc2=scatter!(ax, xy2[1,:],xy2[2,:],color=:red)

Legend(fig[1,2],[sc1, sc2],["nodes on Ω₁","nodes on Ω₂"])

fig
```

## 3D example

The same syntax can be used for regions in three dimensions (and even higher), as the following example shows. 


```@example 3d_sphere
using QuadratureOnImplicitRegions, CairoMakie

ψ(x)= x'*x-1.0 
a,b=zeros(3), ones(3) #the unit cube. 
quad_order=5 

xyz1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)

fig=Figure()
ax = Axis3(fig[1, 1], xlabel = "x", ylabel = "y",zlabel="z",
    title = "Quadrature nodes on Ω₁",
    azimuth = -.2*pi,
    elevation= pi/6 )

scatter!(ax,xyz1[1,:],xyz1[2,:],xyz1[3,:],color=:blue)

fig
```
