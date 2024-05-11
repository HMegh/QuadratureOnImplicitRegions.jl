using LinearAlgebra
using FastGaussQuadrature

sgn=(m,s,σ)-> m==σ*s ? σ*m : 0;

function findroots(ψ_list,domain)
    #Returns a list of union of all zeros of psi_i in  domain and the
    #boundaries of the domain (for convenience). 
    
    a,b=domain;
    Z=[a,b];


    Xsamples=range(a,b,20);
    for ψ in ψ_list
       
        samples=ψ.(Xsamples);

        

        indx=findall(samples[1:end-1] .* samples[2:end] .<0);
        
        zer=
        [bisection(ψ,Xsamples[kk],Xsamples[kk+1])
            for kk in indx];


        Z=vcat(Z,zer);

    end
    return sort(Z);
end
function bisection(f,a::T,b::T) where T<:Union{Float16,Float32,Float64}


    p=f(a);q=f(b);

    #see if the zero is at the endpoints.
    if p==0 return a;end
    if q==0 return b; end;

    if p*q>0 return; end

   

    while abs(b-a)>eps(one(T))
        c=(a+b)/2;r=f(c);



        if r*p<0 
             b=c;q=r;
         else
             a=c;p=r; 
         end
         
        if abs(p)<eps(one(T)) return a end 
        if abs(q)<eps(one(T)) return b end
    end

    return (a+b)/2;
end




function lgwtjl(N::Int,a::T,b::T) where T<:Union{Float16,Float32,Float64}
    xref,wref=gausslegendre(N)
    x=Vector{T}(undef,N)
    w=Vector{T}(undef,N)
    x.=(b-a)/2*xref .+(a+b)/2
    w.=(b-a)/2*wref
    return x,w
end

function linsamples_creator(U::Matrix{T},n::Int) where T<:Union{Float16,Float32,Float64}
    d=size(U,2);

    #linspaces in each direction
    lins=[range(U[1,i],U[2,i],n) for i=1:d];
    
    
    #Iterators.product returns the cartesian product as an object,
    # collect makes it a matrix and [:] makes it a vector of tuples. 
    x= collect(Iterators.product(lins...))[:];


    x=[[y for y in z] for z in x];
    
end   


function tensor_GL_rule(U::Matrix{T},q::Int) where T<:Union{Float16,Float32,Float64}


    d=size(U,2);
    xref, wref=lgwtjl(q,zero(T),one(T));

    X=[U[1,i] .+ xref .* (U[2,i]-U[1,i]) for i=1:d]; 
    W=[ wref .* (U[2,i]-U[1,i]) for i=1:d]; 

    #Iterators.product returns the cartesian product as an object,
    # collect makes it a matrix and [:] makes it a vector of tuples. 
    x=collect(Iterators.product(X...))[:];

    x=[[t for t in y] for y in x]

    w=[prod(A) for A in collect(Iterators.product(W...))[:]];

    return x,w
end
function d1_case(ψ_list,s_list::Vector{T},domain::Array{T},q::Int) where T<:Union{Float16,Float32,Float64}




    #The one-dimensional case. s_list and psi_list
    # have the same length
    zer=findroots(ψ_list,domain);
    x=Vector{T}(undef,0);
    w=Vector{T}(undef,0);

    for i=1:length(zer)-1
        a=zer[i];b=zer[i+1];c=(a+b)/2;
        #check to see if s_i*psi_i(c)>0 for all i      
        #I used count: if there is an instance where 
        # psi_i(c)*s_i<0, then flag is 0. 

        flag=true

        # flag=count([s_list[i]*ψ_list[i](c)<0 
        #         for i=1:length(s_list)])==0;

        for i=1:length(s_list) 
            if s_list[i]*ψ_list[i](c)<0 flag=false end
        end


        if flag #psi_i(c)s_i>0 for all i

            xx,ww=lgwtjl(q,a,b);
            x=vcat(x,xx);
            w=vcat(w,ww);
        end
    end
    

    return x,w
    
end



remove_kth(x::Vector,k::Integer)=vcat(x[1:k-1],x[k+1:end]);
remove_kth(x::Number,k::Integer)=Vector{typeof(x)}(undef,0); 

#tilde insert takes  [1,2,3,4],3,pi and returns [1,2,pi,3,4]
#there is a builtin function insert! but it does not work with scalars
## 
insert_kth(x::Vector,k::Integer,t::Number)=vcat(x[1:k-1],t,x[k:end]);
function insert_kth(x::Number,k::Integer,t::Number)
    if k==1
        return [t,x]
    elseif k==2
        return [x,t]
    end
end





