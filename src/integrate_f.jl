#instead of returning nodes and weights, we can integrate a function f directly.
#this should be faster. 


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
    d=size(U,2)
    n=size(x_ref,1)

    X=[(U[1,i]+U[2,i])/2 .+ x_ref .* (U[2,i]-U[1,i])/2 for i=1:d]
    W=[ w_ref/2 .* (U[2,i]-U[1,i]) for i=1:d]

    x_tuples=collect(Iterators.product(X...))

    W=[prod(A) for A in Iterators.product(W...)]

    return sum(f(x_tuples[i])*W[i] for i in eachindex(x_tuples,W))

end

# d1_int_f(x->x^3,[x->x-1.5],[-1],[0,2],10)≈ 1.5^4/4


function algoim_quad(f::Function,ψ::Function,sgn::Number,a,b,q::Integer)
    
    U=vcat(a',b')

    return algoim_quad(f,[ψ],[∇ψ],[sgn],U,q)
end





function algoim_quad(f,ψ_list,∇ψ_list,s_list,U,q,recursion_depth=1)::Float64

    x_ref,w_ref=gausslegendre(q)
    #integrate a function f 
    d=size(U,2);
    if d==1
        return d1_int_f(f,ψ_list,s_list,U,x_ref,w_ref);
    end


    xc=(U[1,:] + U[2,:])/2

  
    n_samples=100;#less samples: faster

    n_samples=100; #total number of samples 
    nnsamples=Int(round(n_samples^(1/d))) # number of samples per dimension
    Samples=linsamples_creator(U,nnsamples)

    for i=length(s_list):-1:1
        ψ=ψ_list[i];

        ψ_xc=ψ(xc);

        ψ_x=ψ.(Samples) 

        δ=maximum(abs.(ψ_x .- ψ_xc));
        # println("ψ_xc =$ψ_xc  ,δ = $δ")

        #pruning

        # if abs(ψ_xc)>=δ #does not prune certain cases such as norm(x)^2 on [0,1]
        if minimum(ψ_x)*maximum(ψ_x) >=0

            if s_list[i]*ψ_xc>=0
                #remove \psi_i
                s_list=vcat(s_list[1:i-1],s_list[i+1:end]);
                ψ_list=vcat(ψ_list[1:i-1],ψ_list[i+1:end]);
                ∇ψ_list=vcat(∇ψ_list[1:i-1],∇ψ_list[i+1:end]);
            else
                #nothing 
                return 0.0;
            end
        end
    end

    if isempty(s_list)
        x,w=tensor_GL_rule(U,q);
        return w'*f.(x);
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
                    
                    return f(xc)*Volume;

                else
                    return 0.0;
                end
            end

            # println("Domain split into two")

            kk=argmax(U[2,:]-U[1,:]);

            U1=deepcopy(U);U2=deepcopy(U);

            U1[2,kk]=(U[1,kk]+U[2,kk])/2.0;
            U2[1,kk]=(U[1,kk]+U[2,kk])/2.0;
            
            i1=algoim_quad(f,ψ_list,∇ψ_list,s_list,U1,q,recursion_depth+1);
            i2=algoim_quad(f,ψ_list,∇ψ_list,s_list,U2,q,recursion_depth+1);
            return i1+i2; 

        end


        
    end

    U_tilde=hcat(U[:,1:k-1],U[:,k+1:end]);
    
    F_tilde=x->F1d(x,f,ψ_list,s_list,k,xkL,xkU,q)
    return algoim_quad(F_tilde, new_ψ_list,new_∇ψ_list,new_s_list,U_tilde,q)


   

end


function restrict(x::Vector,t::Number,k::Integer)#replace x[k] by t
    return vcat(x[1:k-1],t,x[k+1:end])
end

function restrict(x::Number,t::Number,k::Integer)#replace x[k] by t
   if k==1 return [t,x] end
   if k==2 return [x,t] end
end


function F1d(x,f,ψ_list,s_list,k,a,b,q)

    f_tilde=t->f(restrict(x,t,k))
    ψ_list_tilde=[t->ψ(restrict(x,t,k)) for ψ in ψ_list]

    return d1_int_f(f_tilde,ψ_list_tilde,s_list,[a,b],q)

end

