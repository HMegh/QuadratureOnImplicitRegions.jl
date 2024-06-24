#instead of returning nodes and weights, we can integrate a function f directly.
#this should be faster. 

using RectiGrids #0 allocations cartesian products

"""
    int_f_1d(f,N,a,b,x_ref,w_ref)

    integrate f(x) on [a,b] using a reference quadrature x_ref,w_ref
"""
function int_f_1d(f,a,b,x_ref,w_ref)
    s=0.0;
    for i in eachindex(x_ref,w_ref)
        s+=(b-a)/2*w_ref[i]*f((a+b)/2+x_ref[i]*(b-a)/2)
    end
    s
end

function d1_int_f(f,ψ_list,s_list,domain,x_ref,w_ref)
    #The one-dimensional case. s_list and psi_list
    # have the same length
    Z=find_roots(ψ_list,domain[1],domain[2])
    
    s=0.0

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
            s+=int_f_1d(f,a,b,x_ref,w_ref)
        end
    end 

    return s
    
end
"""
    tensor_GL_int_f(f,U,x_ref,w_ref)

Integrate a function f on a cuboid U using a reference quadrature x_ref,w_ref
"""

function tensor_GL_int_f(f,U,x_ref,w_ref)
  
  
    
    X=[muladd.(x_ref,(U[2,i]-U[1,i])/2,(U[1,i]+U[2,i])/2) for i in axes(U,2)]

    W= [ (U[2,i]-U[1,i])/2 *w_ref for i in axes(U,2)]

    x_tuples=grid(X...)
    w_tuples=grid(W...)


    return sum(f(x_tuples[i])*prod(w_tuples[i]) for i in eachindex(x_tuples,w_tuples))

end

# d1_int_f(x->x^3,[x->x-1.5],[-1],[0,2],10)≈ 1.5^4/4


function algoim_quad(f,ψ,sgn,a,b,q)
    U=vcat(a',b')
    x_ref,w_ref=gausslegendre(q)


    return algoim_quad(f,[ψ],[sgn],U,x_ref,w_ref)
end





function algoim_quad(f::F,ψ_list,s_list::Vector{T},U::M,x_ref,w_ref,recursion_depth=1) where {T<:Number,M<:AbstractMatrix,F<:Function}

    #integrate a function f 
    d=size(U,2);
    if d==1
        return d1_int_f(f,ψ_list,s_list,U,x_ref,w_ref)
    end

    xc=(U[1,:] + U[2,:])/2 #midpoint

    n_samples=100; #total number of samples 
    nnsamples=Int(round(n_samples^(1/d))) # number of samples per dimension
    Samples=linsamples_creator(U,nnsamples)

    for i=length(s_list):-1:1
        ψ=ψ_list[i]
        min_ψ=minimum(ψ(x) for x in eachcol(Samples))
        max_ψ=maximum(ψ(x) for x in eachcol(Samples))

        #= 
        Here, I deviated a bit from the paper, to make the case ψ(x)=norm(x)^2 works on [0,ε]^d 
        =#
        if min_ψ*max_ψ≥0
            if s_list[i]*ψ(xc)≥0 
                s_list=deleteat!(copy(s_list),i)
                ψ_list=deleteat!(copy(ψ_list),i)
            else
                return 0.0#TODO: change 0.0 to zero(T)
            end
        end

    end

    if isempty(s_list)
        return tensor_GL_int_f(f,U,x_ref,w_ref)
    end
        
    new_ψ_list=Vector{Function}(undef,0)
    new_s_list=Vector{Float64}(undef,0)

    #the direction of the most change in ψ_1
    ψ_1=ψ_list[1];
    k=argmax(abs.(ForwardDiff.gradient(ψ_1,xc)))

 
    xkL=U[1,k];
    xkU=U[2,k];

    g=zeros(d)
    gtemp=zeros(d)
    δ=fill(-Inf,d)
    for i in eachindex(s_list)
        ψ=ψ_list[i]

        # δ[:] .=fill(-Inf,d)
        fill!(δ,-Inf)

        # g[:] .=ForwardDiff.gradient(ψ,xc)
        finite_difference_gradient!(g,ψ,xc)

        for i in axes(Samples,2)
            # gtemp[:] .=ForwardDiff.gradient(ψ,Samples[:,i])
            finite_difference_gradient!(gtemp,ψ, @view Samples[:,i])
            for kk=1:d 
                δ[kk]=max(δ[kk],abs(gtemp[kk]-g[kk]))
            end
        end
    

        if abs(g[k])>δ[k]&& norm(g + δ)^2/(g[k]-δ[k])^2<20

            ψ_L=(x->ψ(insert_kth(x,k,xkL)))
            ψ_U=(x->ψ(insert_kth(x,k,xkU)))

            si_L=sgn(g[k],s_list[i],-1);
            si_U=sgn(g[k],s_list[i],1);


            # new_ψ_list=vcat(new_ψ_list,ψ_L,ψ_U); 
            # new_s_list=vcat(new_s_list,si_L,si_U)
            push!(new_ψ_list,ψ_L,ψ_U)
            push!(new_s_list,si_L,si_U)

        else #subdivide the domain (unless it is too small)
            Volume=prod(U[2,:] - U[1,:]);

            if recursion_depth>=16
                flag=true
                for i in eachindex(s_list)
                    if s_list[i]*ψ_list[i](xc)<= 0 flag=false end
                end

                printstyled("Warning: lower order method used in $U\n";color=:red)

                if flag
                    return f(xc)*Volume 
                else
                    return 0.0
                end
            end 

            kk=argmax(U[2,:]-U[1,:])
            U1=copy(U);U2=copy(U);

            U1[2,kk]=(U[1,kk]+U[2,kk])/2.0;
            U2[1,kk]=(U[1,kk]+U[2,kk])/2.0;



            I1=algoim_quad(f,ψ_list,s_list,U1,x_ref,w_ref,recursion_depth+1)
            I2=algoim_quad(f,ψ_list,s_list,U2,x_ref,w_ref,recursion_depth+1)

          
            
           
            return I1+I2
        end

    end


    U_tilde=hcat(U[:,1:k-1],U[:,k+1:end]);
    
    F_tilde(x)= F1d(x,f,ψ_list,s_list,k,xkL,xkU,x_ref,w_ref)
    return algoim_quad(F_tilde, new_ψ_list,new_s_list,U_tilde,x_ref,w_ref)


   

end


function restrict(x::V,t::N,k::I) where {V<:AbstractVector,N,I}#replace x[k] by t
    # return vcat(x[1:k-1],t,x[k+1:end])
    z=copy(x)
    z[k]=t
    return z
end

"""
    restrict!(x::V,t::N,k::I,z_temp::V) where {V<:AbstractVector,N,I}
    x: vector of size d.
    z: vector of size d. 

    z[1:k-1] =x[1:k-1]
    z[k]=t 
    z[k+1:end]= x[k:end]
"""

function restrict!(z,x::V,t::N,k::I) where {V<:AbstractVector,N,I}
    for i=1:k-1 z[i]=x[i] end 
    z[k]=t
    for i=k:length(x) z[i+1]=x[i] end
    return z
end

function restrict(x::V,t::N,k::I) where {V<:Number,N,I}#replace x[k] by t
    # return vcat(x[1:k-1],t,x[k+1:end])
   if k==1 return [t;x] end 
   return [x;t] 
end

function restrict!(z::V,x::W,t::N,k::I) where {V<:AbstractVector,N,I,W<:AbstractVector}
    if k==1 
        z[1]=t;z[2]=x;
    else
        z[1]=x;z[2]=t;
    end
    return z
end



function F1d(x,f,ψ_list,s_list,k,a,b,x_ref,w_ref)

    f_tilde=t->f(restrict(x,t,k))
    ψ_list_tilde=[t->ψ(restrict(x,t,k)) for ψ in ψ_list]

    return d1_int_f(f_tilde,ψ_list_tilde,s_list,[a,b],x_ref,w_ref)

end

