

function  d1_case_only_bounds(ψ_list,s_list,domain,order)
    Z=find_roots(ψ_list,domain[1],domain[2])
    
    #estimate how many b
    nbr_subint= d1_count_subintervals(ψ_list,s_list,domain,Z)

    V=Vector{Tuple{Float64,Float64}}(undef,nbr_subint)


    n=length(Z)
    idx=1
    for i=1:n-1
        a=Z[i];b=Z[i+1];c=(a+b)/2

        flag=true 
        for i in eachindex(ψ_list) 
            if s_list[i]*ψ_list[i](c)<0
                flag=false;break;
            end
        end    
        if flag
            V[idx]= (a,b)
            idx+=1
        end

    end 

    return Imp_d_1_rule(V,order)
end