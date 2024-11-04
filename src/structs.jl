
"""
    struct Imp_d_1_rule{T} 

Stores the sub-intervals on which the quadrature points are defined.
Fields: 
- `subintervals::Vector{Tuple{T,T}}`
- `order::Int`

example: if ``ψ(x)=x^2-2``, ``Ω=[-2,2]`` and we are looking for the subset ``ψ(x)>0`, then the `subintervals` are ``(-2,-√2), (√2,2)``. 

The order is the order of quadrature rule used. 
"""

import Base:collect

struct Imp_d_1_rule{T} 
    subintervals::Vector{Tuple{T,T}}
    order::Int
end


Q=Imp_d_1_rule([(-2.0,-√2),(√2,2.0)],10)

function collect(Q::Imp_d_1_rule{T}) where T
    (; subintervals,order) =Q 
    xg,wg=gausslegendre(order)

    nodes  =Vector{T}(undef,length(subintervals)*order)
    weights=Vector{T}(undef,length(subintervals)*order)
    for j in eachindex(subintervals)
        a,b=subintervals[j]
        for k=1:order
            nodes[(j-1)*order+k]  =xg[k]*(b-a)/2+(a+b)/2
            weights[(j-1)*order+k]=wg[k]*(b-a)/2
        end 
    end
    return nodes,weights

end


function collect(Q::Imp_d_1_rule{T},xg::Vector{T},wg::Vector{T}) where T
    (; subintervals,order) =Q 
    @assert length(xg)== length(wg)==order 

    nodes  =Vector{T}(undef,length(subintervals)*order)
    weights=Vector{T}(undef,length(subintervals)*order)
    for j in eachindex(subintervals)
        a,b=subintervals[j]
        for k=1:order
            nodes[(j-1)*order+k]  =xg[k]*(b-a)/2+(a+b)/2
            weights[(j-1)*order+k]=wg[k]*(b-a)/2
        end 
    end
    return nodes,weights

end


"""


x,w= collect(Q)

sum(w) -(4-2sqrt(2))



xg,wg=gausslegendre(10)

@btime collect(Q,xg,wg)
#   47.891 ns (5 allocations: 480 bytes)
"""