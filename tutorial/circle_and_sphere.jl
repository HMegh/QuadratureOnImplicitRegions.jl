using QuadratureOnImplicitRegions
using Test
using Plots #to plot the quad points


#Example 1: 
#consider the unit interval Ω=[0,1]², split into two regions:
#Ω₁= {(x,y)∈Ω | x^2+y^2< 1 }
# and Ω₂= {(x,y)∈Ω | x^2+y^2>1 }
# Our goal is to generate quadrature nodes and weights on each region (subdomain)

ψ(x)= x'*x-1.0 #basically x^2+y^2 -1 
a,b=zeros(2), ones(2) #the unit interval. 
quad_order=10 #the higher, the more accurate the quadrature is

xy1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)
xy2,w2=algoim_nodes_weights(ψ,+1.0, a,b,quad_order)

#test the sum of the weights (it is equal to the area of each region)
@test sum(w1)≈ π/4 
@test sum(w2)≈ (1-π/4)

#xy a vector of vectors:
x1,y1= first.(xy1) , last.(xy1) 
x2,y2= first.(xy2) , last.(xy2) 

#visualize the quadrature nodes:

#the rectangle
rect_x=[a[1],b[1],b[1],a[1],a[1]]
rect_y=[a[2],a[2],b[2],b[2],a[2]]
P=plot(rect_x,rect_y,c=:black,label="Ω",legend=:outerleft,axis=false,grid=false)

#the interface ψ=0
plot!(P,range(a[1],b[1],50),
range(a[2],b[2],50),(x,y)->ψ([x,y]),levels=[0.0],c=:green,seriestype=:contour,label="ψ=0",cbar=false)

#the quadrature nodes
scatter!(P,x1,y1,c=:blue,label="Ω₁")
scatter!(P,x2,y2,c=:red,label="Ω₂")
savefig(P,"tutorial/example_1.png")



#Example 2: 
#Consider the domain [-2,2]² split by the unit circle into the unit disk and the remainder. 
#It can be handled similarly to the previous example. The algorithm will detect the need to split the domain automatically. 


ψ(x)= x'*x-1.0 #the same as before 
a,b=-2ones(2), 2ones(2) 


quad_order=5  


xy1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)
xy2,w2=algoim_nodes_weights(ψ,+1.0, a,b,quad_order)

# need quad_order≥20 to achieve this accuracy. For the sake of this tutorial, we will use a smaller quad_order
# @test sum(w1)≈ π
# @test sum(w2)≈ 16-π


#plotting the quad points as before
x1,y1= first.(xy1) , last.(xy1) 
x2,y2= first.(xy2) , last.(xy2) 


rect_x=[a[1],b[1],b[1],a[1],a[1]]
rect_y=[a[2],a[2],b[2],b[2],a[2]]
P2=plot(rect_x,rect_y,c=:black,label="Ω",legend=:outerleft,axis=false,grid=false)

plot!(P2,range(a[1],b[1],50),
range(a[2],b[2],50),(x,y)->ψ([x,y]),levels=[0.0],c=:green,seriestype=:contour,label="ψ=0",cbar=false)
scatter!(P2,x1,y1,c=:blue,label="Ω₁")
scatter!(P2,x2,y2,c=:red,label="Ω₂")
savefig(P2,"tutorial/example_2.png")


#Example 3:
#It is the same as example 1, but in three dimensions



ψ(x)= x'*x-1.0 #the same as before 
a,b=zeros(3), ones(3) 

quad_order=5  

xyz1,w1=algoim_nodes_weights(ψ,-1.0, a,b,quad_order)
xyz2,w2=algoim_nodes_weights(ψ,+1.0, a,b,quad_order)

x1,y1,z1= first.(xyz1),map(t->t[2],xyz1) , last.(xyz1) 
x2,y2,z2= first.(xyz2) ,map(t->t[2],xyz2),  last.(xyz2) 

#the unit cube
A=[0,1,1,0,0]
B=[0,0,1,1,0]
C=[0,0,0,0,0]
D=[1,1,1,1,1]
P3=plot()
plot!(P3,A,B,C,label="Ω",c=:black,axis=false,grid=false)
plot!(P3,A,B,D,c=:black,label=nothing)
plot!(P3,A,C,B,c=:black,label=nothing)
plot!(P3,A,D,B,c=:black,label=nothing)

P31=deepcopy(P3)
P32=deepcopy(P3)
scatter!(P31,x1,y1,z1,c=:blue,label="Ω₁")
scatter!(P32,x2,y2,z2,c=:red,label="Ω₂")



plot!(P31,camera=(60,15))

# plot(P32,camera=(30,30))

savefig(P31,"tutorial/example_3.png")
