


include("algoim_functions.jl")


import ForwardDiff


function algoim_nodes_weights(ψ::Function,sgn::Number,a::Vector{Float64},b::Vector{Float64},q::Int) 

    U=vcat(a',b')
    ∇ψ = x-> ForwardDiff.gradient(ψ,x)
    return algoim_nodes_weights([ψ],[∇ψ], [sgn],U,q)
end




function algoim_nodes_weights(ψ_list::Vector,∇ψ_list::Vector,s_list::Vector,U,q::Integer,recursion_depth=1)


    d=size(U,2);
    if d==1
        x,w=d1_case(ψ_list,s_list,U,q);
        return x,w;
    end


    xc=(U[1,:] + U[2,:])/2

    n_samples=15625;
    n_samples=4096;#less samples: faster
    n_samples=729;#less samples: faster

    nnsamples=Int(round(n_samples^(1/d)));

    Samples=linsamples_creator(U,nnsamples)

    for i=length(s_list):-1:1
        ψ=ψ_list[i];

        ψ_xc=ψ(xc);

        ψ_x=ψ.(Samples) 

        δ=maximum(abs.(ψ_x .- ψ_xc));
        # println("ψ_xc =$ψ_xc  ,δ = $δ")

        #pruning

        # if abs(ψ_xc)>=δ #does not prune certain cases such that norm(x)^2 on [0,1]
        if minimum(ψ_x)*maximum(ψ_x) >=0

            if s_list[i]*ψ_xc>=0
                #remove \psi_i
                s_list=vcat(s_list[1:i-1],s_list[i+1:end]);
                ψ_list=vcat(ψ_list[1:i-1],ψ_list[i+1:end]);
                ∇ψ_list=vcat(∇ψ_list[1:i-1],∇ψ_list[i+1:end]);
            else
                #nothing 
                return Vector{Float64}[],Float64[];
            end
        end
    end

    if isempty(s_list)
        x,w=tensor_GL_rule(U,q);
        return x,w;
    end


    new_ψ_list=Vector{Function}(undef,0)
    new_∇ψ_list=Vector{Function}(undef,0)
    new_s_list=Vector{Float64}(undef,0)

        

    ψ_1=ψ_list[1];
    ∇ψ_1=∇ψ_list[1];

    #this does not strike me as the best method to choose k
    k=argmax(abs.(∇ψ_1(xc)))
    # println("dimension $(length(xc)), direction chosen $k, U=$U")
    
    # k=argmax(abs.(sum(∇ψ_1.(Samples))))

    xkL=U[1,k];
    xkU=U[2,k];

    for i=1:length(s_list)
        ψ=ψ_list[i]
        ∇ψ=∇ψ_list[i]

        g=∇ψ(xc);

        ∇ψ_x=hcat([∇ψ(Samples[ii])  for ii=1:n_samples]...)

        
        #\delta is the max on each row
        δ=[maximum(abs.(∇ψ_x[kk,:] .- g[kk])) for kk=1:d]


        if abs(g[k])>δ[k]&& norm(g + δ)^2/(g[k]-δ[k])^2<20



            #The restriction to the faces in the e_k direction.

            #I am using  [x][1:k-1] instead of x[1:k-1] since 
            # Vector[1:0] is well defined (empty array) but Number[1:0] is not. 
            

            ψ_L=(x->ψ(insert_kth(x,k,xkL)))
            ψ_U=(x->ψ(insert_kth(x,k,xkU)))

            ∇ψ_L=x->remove_kth(∇ψ(insert_kth(x,k,xkL)),k)
            ∇ψ_U=x->remove_kth(∇ψ(insert_kth(x,k,xkU)),k)
            
                

            # println("xkL,xkU=$xkL, $xkU , psi L and psi U defined")

            

            si_L=sgn(g[k],s_list[i],-1);
            si_U=sgn(g[k],s_list[i],1);

            new_ψ_list=vcat(new_ψ_list,ψ_L,ψ_U); 
            new_∇ψ_list=vcat(new_∇ψ_list,∇ψ_L,∇ψ_U); 
            new_s_list=vcat(new_s_list,si_L,si_U)

        else

            Volume=prod(U[2,:] - U[1,:]);

            # if Volume<1e-3 #should be changed to something like #iterations>15
              if recursion_depth>=16

                #flag=1 if psi_i*s_i>0 for all i

                flag=prod([s_list[i]*ψ_list[i](xc)>0 for i=1:length(s_list)]);
                # display(flag)
                # display(∇ψ_x)
                # println("xc= $xc, ψ=[ $(ψ_list[1](xc))]")
                # println(" ψ(0,0)= $(ψ_list[1]([0,0.0]))")
                 printstyled("Warning: lower order method used in $U\n";color=:red)

                #  println((g,δ,k,abs(g[k]), δ[k], norm(g + δ)^2/(g[k]-δ[k])^2))
                if flag #low order method
                    
                    return xc, Volume;

                else
                    return Vector{Float64}[],Float64[];
                end
            end

            # println("Domain split into two")

            kk=argmax(U[2,:]-U[1,:]);

            U1=deepcopy(U);U2=deepcopy(U);

            U1[2,kk]=(U[1,kk]+U[2,kk])/2.0;
            U2[1,kk]=(U[1,kk]+U[2,kk])/2.0;
            
            x1,w1=algoim_nodes_weights(ψ_list,∇ψ_list,s_list,U1,q,recursion_depth+1);
            x2,w2=algoim_nodes_weights(ψ_list,∇ψ_list,s_list,U2,q,recursion_depth+1);

            x=vcat(x1,x2); w=vcat(w1,w2);

            return x,w; 

        end


        
    end

    U_tilde=hcat(U[:,1:k-1],U[:,k+1:end]);

    x_tilde,w_tilde=algoim_nodes_weights(new_ψ_list,new_∇ψ_list,new_s_list,U_tilde,q);

    x=Vector{Float64}[];
    w=Float64[];

    for ii=1:length(x_tilde)

        QQ=x_tilde[ii] # a (d-1) tuple
        PP=w_tilde[ii]
        ψ_list_tilde=
        [t->ψ(insert_kth(QQ,k,t)) for ψ in ψ_list];
        
        x_slice,w_slice=d1_case(ψ_list_tilde,s_list,[xkL;xkU],q);

        

        if ~isempty(x_slice)
            
            w=vcat(w,PP*w_slice)


            # println((PP,sum(w_slice),length(w),sum(w)))



            #QQ has d-1 components, we want to insert x_slice in the k-th components


            x_loc=[insert_kth(QQ,k,z) for z in x_slice]
            
            x=vcat(x,x_loc);

        end

    end





    return x,w;


end


