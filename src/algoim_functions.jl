using LinearAlgebra
using FastGaussQuadrature

sgn=(m,s,σ)-> m==σ*s ? σ*m : 0;

function findroots(ψ_list,domain)
    #Returns a list of union of all zeros of psi_i in  domain and the
    #boundaries of the domain (for convenience). 
    
    a,b=domain;
    Z=[a,b];


    Xsamples=range(a,b,20);
    for jj=1:length(ψ_list)
        ψ=ψ_list[jj];
       
        samples=ψ.(Xsamples);

        

        indx=findall(samples[1:end-1] .* samples[2:end] .<0);
        
        zer=
        [bisection(ψ,Xsamples[kk],Xsamples[kk+1])
            for kk in indx];


        Z=vcat(Z,zer);

    end
    return sort(Z);
end
function bisection(f,a,b)


    p=f(a);q=f(b);

    #see if the zero is at the endpoints.
    if p==0 return a;end
    if q==0 return b; end;

    if p*q>0 return; end

   

    while abs(b-a)>eps(1.0)
        c=(a+b)/2;r=f(c);



        if r*p<0 
             b=c;q=r;
         else
             a=c;p=r; 
         end
         
        if abs(p)<eps(1.0) return a end 
        if abs(q)<eps(1.0) return b end
    end

    return (a+b)/2.0;
end


#The following functions is much slower than the one in FastGaussQuadrature
# function lgwtjl(N,a,b)

#     #adaptation of lgwt for MATLAB

#     N1=N+1; N2=N+2;

#     xu=range(-1,1,N1);

#     y=cos.((2*(0:N) .+ 1)*pi/(2 * N+2))+(0.27/N1).*sin.(pi*xu*N/N2);


#     L=zeros(N1,N2);
#     Lp=zeros(N1,N2);


#     y0=ones(N1);
#     Lq=y0;

#     while maximum(abs.(y .- y0))>eps(Float64)
        
        
#         L[:,1].=1;
#         Lp[:,1].=0;
        
#         L[:,2]=y;
#         Lp[:,2].=1;
        
#         for k=2:N1
#             L[:,k+1]=( (2*k-1)*y.*L[:,k]-(k-1)*L[:,k-1] )/k;
#         end
    
#         Lq=(N2)*( L[:,N1]-y.*L[:,N2] )./(1 .- y.^2);   
#         y0=y;
#         y=y0 -L[:,N2]./Lq;
        
#     end

#     x=(a*(1.0 .-y).+ b*(1.0 .+y))/2;      

#     w=(b-a)./((1 .- y.^2).*Lq.^2)*(N2/N1)^2;
#     return x,w;

# end



function lgwtjl(N,a,b)
    x,w=gausslegendre(N)
    x=(b-a)/2*x .+(a+b)/2
    w=(b-a)/2*w
    return x,w
end

function linsamples_creator(U,n)
    d=size(U,2);

    #linspaces in each direction
    lins=[range(U[1,i],U[2,i],n) for i=1:d];
    
    
    #Iterators.product returns the cartesian product as an object,
    # collect makes it a matrix and [:] makes it a vector of tuples. 
    x= collect(Iterators.product(lins...))[:];


    x=[[y for y in z] for z in x];
    
end   


function tensor_GL_rule(U,q)


    d=size(U,2);
    xref, wref=lgwtjl(q,0,1);

    X=[U[1,i] .+ xref .* (U[2,i]-U[1,i]) for i=1:d]; 
    W=[ wref .* (U[2,i]-U[1,i]) for i=1:d]; 

    #Iterators.product returns the cartesian product as an object,
    # collect makes it a matrix and [:] makes it a vector of tuples. 
    x=collect(Iterators.product(X...))[:];

    x=[[t for t in y] for y in x]

    w=[prod(A) for A in collect(Iterators.product(W...))[:]];

    return x,w
end
function d1_case(ψ_list,s_list,domain,q)




    #The one-dimensional case. s_list and psi_list
    # have the same length
    zer=findroots(ψ_list,domain);
    x=Vector{Float64}(undef,0);
    w=Vector{Float64}(undef,0);
    flag=1;

    for i=1:length(zer)-1
        a=zer[i];b=zer[i+1];c=(a+b)/2;
        #check to see if s_i*psi_i(c)>0 for all i      
        #I used count: if there is an instance where 
        # psi_i(c)*s_i<0, then flag is 0. 

        flag=count([s_list[i]*ψ_list[i](c)<0 
                for i=1:length(s_list)])==0;

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





