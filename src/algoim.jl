


function algoim_nodes_weights(ψ::F,sgn::N,a::Vector{T},b::Vector{S},q::I) where {F,N,T,I,S}
    x_ref,w_ref=gausslegendre(q)
    U=Float64.(vcat(a',b'))
    return algoim_nodes_weights([ψ], [sgn],U,x_ref,w_ref)
end

function algoim_nodes_weights(ψ_list::Vector{F},s_list::Vector{S},U::Matrix{T},q::I) where {T,F,S,I<:Integer}
    x_ref,w_ref=gausslegendre(q)
    U=Float64.(vcat(a',b'))
    return algoim_nodes_weights(ψ_list,s_list,U,x_ref,w_ref)
end

function algoim_nodes_weights(ψ_list,s_list,U::Matrix{T},x_ref,w_ref,recursion_depth=1) where T

    


    d=size(U,2);
    if d==1
        x,w=d1_case(ψ_list,s_list,U,x_ref,w_ref)
        return x,w
    end

    xc=(U[1,:] + U[2,:])/2 #midpoint

    n_samples=100; #total number of samples 
    nnsamples=Int(round(n_samples^(1/d))) # number of samples per dimension
    Samples=linsamples_creator(U,nnsamples)

    for i=length(s_list):-1:1
        ψ=ψ_list[i]
        ψ_xc=ψ(xc)
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
                return Matrix{T}(undef,d,0), Float64[]
            end
        end

    end

    if isempty(s_list)
        x,w=tensor_GL_rule(U,x_ref,w_ref)
        return x,w;
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
                    return xc, Volume 
                else
                    return Vector{Float64}[],Float64[];
                end
            end 

            kk=argmax(U[2,:]-U[1,:])
            U1=copy(U);U2=copy(U);

            U1[2,kk]=(U[1,kk]+U[2,kk])/2.0;
            U2[1,kk]=(U[1,kk]+U[2,kk])/2.0;



            x1,w1=algoim_nodes_weights(ψ_list,s_list,U1,x_ref,w_ref,recursion_depth+1);
            x2,w2=algoim_nodes_weights(ψ_list,s_list,U2,x_ref,w_ref,recursion_depth+1);

          
            
            x= hcat(x1,x2)
            w=vcat(w1,w2)
            return x,w
        end

    end

    U_tilde=hcat(U[:,1:k-1],U[:,k+1:end])
    
    x_tilde,w_tilde=algoim_nodes_weights(new_ψ_list,new_s_list,U_tilde,x_ref,w_ref) #(d-1)xN  matrix 


    QQ=Vector{T}(undef,d-1)

    #This part is the one using 90% of the allocations 

    #estimate the number of nodes
    cnt=0
    for ii in axes(x_tilde,2)
        QQ[:] .=@view x_tilde[:,ii] #(d-1) vector
        # PP=w_tilde[ii]
        ψ_list_tilde=[t->ψ(insert_kth(QQ,k,t)) for ψ in ψ_list]
        cnt+=d1_count_subintervals(ψ_list_tilde,s_list,[xkL;xkU])


      
    end

    x=Matrix{T}(undef,d,cnt*length(x_ref))
    w=Vector{T}(undef,cnt*length(x_ref))


    j=1
    for ii in axes(x_tilde,2)
        QQ[:] .=@view x_tilde[:,ii]
        PP=w_tilde[ii]

        ψ_list_tilde=[t->ψ(insert_kth(QQ,k,t)) for ψ in ψ_list]
        x_slice,w_slice=d1_case(ψ_list_tilde,s_list,[xkL;xkU],x_ref,w_ref)
        
        n=length(w_slice)

        if isempty(w_slice) continue end

        w[j:j+n-1] .= PP*w_slice 
        # x[:,j:j+n] = [repeat(QQ[1:k-1],1,n);x_slice;repeat(QQ[k:end],1,n)]
        x[1:k-1,j:j+n-1] = repeat(QQ[1:k-1],1,n)
        x[k,j:j+n-1] = x_slice
        x[k+1:end,j:j+n-1] = repeat(QQ[k:end],1,n)
        j+=n

    end


    return x,w

end