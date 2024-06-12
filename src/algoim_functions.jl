sgn(m,s,σ)= ( m==σ*s ? σ*m : 0.0)


find_root(ψ::F,a::T,b::T) where {F,T<:Integer} = find_root(ψ,Float64(a),Float64(b))


"""
    find_root(ψ::F,a::T,b::T) 
    
Finds a root of a function ψ in the open interval `(a,b)`
"""
function find_root(ψ::F,a::T,b::T) where {F,T<:AbstractFloat}
    if ψ(a)*ψ(b)≥0 error("Find_root failed: Choose a smaller interval") end 

    c=b
    cnext=(a+b)/2
    while abs(c-cnext)>4eps(T) 
        c=cnext
        cnext=cnext-ψ(cnext)/ForwardDiff.derivative(ψ,cnext)
    end
    return c
end

"""
    find_roots(ψ::F,a::T,b::T) where {F,T}

Returns multiple roots of ψ in the open interval `[a,b]`
"""
function find_roots(ψ::F,a::T,b::T) where {F,T}

    n=20 #split the interval [a,b] into n sub-intervals
    h=(b-a)/n
    Z=Vector{Float64}(undef,0)
    for i=0:n-1
        if ψ(a+i*h)*ψ(a+(i+1)*h)<0 
            z=find_root(ψ,a+i*h,a+(i+1)*h)
            push!(Z,z)
        end
    end
    
    for i=1:n if ψ(a+i*h)==0 push!(Z,a+i*h) end end 

    sort!(Z)

    #remove duplicates
    for i=length(Z):-1:2
        if abs(Z[i]-Z[i-1])≤4eps(T)
            deleteat!(Z,i)
        end
    end

    return Z

end


"""
    find_roots(ψ_list::Vector{F},a::T,b::T) where {F,T}

Returns multiple roots of multiple functions ψ in the closed interval `[a,b]`, it includes the endpoints by default even if they are not roots. 
"""
function find_roots(ψ_list::Vector{F},a::T,b::T) where {F,T}

    n=20 #split the interval [a,b] into n sub-intervals
    h=(b-a)/n
    Z=[a,b]
    for ψ in ψ_list
        zer=find_roots(ψ,a,b)
        append!(Z,zer)
    end
    sort!(Z)
    Z
end

"""
    shifted_gl(x_ref::Vector{T},w_ref::Vector{T},a::T,b::T)

Shifts the Gauss-Legendre rule (x_ref,w_ref) from `(-1,1)` to the interval [a,b]
"""
function shifted_gl(x_ref::Vector{T},w_ref::Vector{T},a::T,b::T) where T
    x,w=similar(x_ref),similar(w_ref)
    for i in eachindex(x) x[i]=(a+b)/2+(b-a)/2*x_ref[i] end
    for i in eachindex(w) w[i]=(b-a)/2*w_ref[i] end
    return x,w
end

"""
    shifted_gl!(x_ref::Vector{T},w_ref::Vector{T},a::T,b::T,x::Vector{T},w::Vector{T})

Shifts the Gauss-Legendre rule (x_ref,w_ref) from `(-1,1)` to the interval [a,b] and stores them in x,w.
"""
function shifted_gl!(x_ref::Vector{T},w_ref::Vector{T},a::T,b::T,x::Vector{T},w::Vector{T}) where T
    for i in eachindex(x) x[i]=(a+b)/2+(b-a)/2*x_ref[i] end
    for i in eachindex(w) w[i]=(b-a)/2*w_ref[i] end
    return 
end


"""
    d1_case(ψ_list,s_list,domain,x_ref,w_ref)

    Returns a quadrature rule on the subset where s_i*ψ_i<0 (see the paper).
"""
function  d1_case(ψ_list,s_list,domain,x_ref,w_ref)
    Z=find_roots(ψ_list,domain[1],domain[2])
    x=Vector{eltype(x_ref)}(undef,0)
    w=Vector{eltype(w_ref)}(undef,0)

    xloc=similar(x_ref)
    wloc=similar(w_ref)

    n=length(Z)
    for i=1:n-1
        a=Z[i];b=Z[i+1];c=(a+b)/2

        flag=true 
        for i in eachindex(ψ_list) 
            if s_list[i]*ψ_list[i](c)<0
                flag=false;break;
            end
        end    
        if flag
            shifted_gl!(x_ref,w_ref,a,b,xloc,wloc)
            append!(x,xloc)
            append!(w,wloc)
        end
    end 

    return reshape(x,1,:),w
end


"""
    tensor_GL_rule(U::Matrix{T},x_ref::Vector{T},w_ref::Vector{T})

    Returns a matrix x and a vector w (multidimensional quadrature) 

    Credit for the idea: 
    https://discourse.julialang.org/t/about-storing-the-results-of-iterators-product-into-an-array-efficiently/115504
"""
function tensor_GL_rule(U::Matrix{T},x_ref::Vector{T},w_ref::Vector{T}) where T
    d=size(U,2)
    n=size(x_ref,1)

    X=[(U[1,i]+U[2,i])/2 .+ x_ref .* (U[2,i]-U[1,i])/2 for i=1:d]
    W=[ w_ref/2 .* (U[2,i]-U[1,i]) for i=1:d]

    x_tuples=collect(Iterators.product(X...))

    x=Matrix(reshape(reinterpret(T, x_tuples), d, n^d))

    
    w=[prod(A) for A in collect(Iterators.product(W...))[:]]

    return x,w
end

"""
    linsamples_creator(U::Matrix{T},n::I) where {T,I}
    
    returns sample points in the cuboid U (as columns).
"""

function linsamples_creator(U::Matrix{T},n::I) where {T,I}
    d=size(U,2);

    #linspaces in each direction
    lins=[range(U[1,i],U[2,i],n) for i=1:d];
    
    
    #Iterators.product returns the cartesian product as an object,
    # collect makes it a matrix and [:] makes it a vector of tuples. 
    x_tuples= collect(Iterators.product(lins...))


    return Matrix(reshape(reinterpret(T, x_tuples), d, n^d))
    
end  


# remove_kth(x::Vector,k::Integer)=vcat(x[1:k-1],x[k+1:end]);
# remove_kth(x::Number,k::Integer)=Vector{typeof(x)}(undef,0); 

function remove_kth(x::V,k::I) where {I,V<:AbstractVector}
    return deleteat!(copy(x),k)
end

function remove_kth(x::T,k::I) where {I,T<:Number}
    return Vector{T}(undef,0)
end

function insert_kth(x::V,k::I,t::T) where {I,V<:AbstractVector,T<:Number}
    [x[1:k-1];t;x[k:end]]
end

function insert_kth(x::N,k::I,t::S) where {N<:Number,I,S<:Number}
    if k==1
        return [t,x]
    elseif k==2
        return [x,t]
    end
end
