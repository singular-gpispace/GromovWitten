@doc raw"""
    feynman_integral_branch_type(  F::FeynmanGraph,   a::Vector{Int64} ;l=zeros(Int,nv(F.G)))

 compute the Feynman Integral for a specified branch type `a` for all ordering `Ω`
    
# Examples (without vertex contribution)
 
```julia
julia> G=FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> F=FeynmanIntegral(G);

julia> a=[0,2,1,0,0,1];

julia> feynman_integral_branch_type(F,a)
256*q[2]^4*q[3]^2*q[6]^2
```

     feynman_integral_branch_type( F::FeynmanGraph,   a::Vector{Int64} ;aa=0,l=zeros(Int,nv(F.G)),g=zeros(Int,nv(F.G)))

# Examples (with vertex contribution)

```julia
julia> G=FeynmanGraph([(1, 2), (2, 3), (1, 3)])
FeynmanGraph([(1, 2), (2, 3), (1, 3)])

julia> F=FeynmanIntegral(G);

julia> a=[0,0,3];
julia> g=[1,0,0];

julia> feynman_integral_branch_type(F,a,g)
115//3*q[3]^6
```
"""
function feynman_integral_branch_type(F::FeynmanIntegral, a::Vector{Int64}; l=zeros(Int, nv(F.G)))
    ee = Edge.(F.G.edge)
    N = sum(a)
    f = signature_and_multiplicities(F.G, a)
    p = 0
    L = lis(F.G, N, l)
    x = F.S[2]  # Vector of elements x from the polynomial ring R
    q = F.S[3]  # Vector of elements q from the polynomial ring R

    for i in eachindex(f)
        tmp = 1  # Use the identity element of the polynomial ring R

        for j in eachindex(f[i][2])
            if f[i][2][j] == -1
                tmp *= constterm(x[src(ee[j])], x[dst(ee[j])], N)
            elseif f[i][2][j] == 0
                tmp *= constterm(x[dst(ee[j])], x[src(ee[j])], N)
            elseif f[i][2][j] == -2
                tmp *= loopterm(q[j], a[j])
            else
                tmp *= proterm(x[src(ee[j])], x[dst(ee[j])], q[j], f[i][2][j], N)
            end
        end
        p += f[i][1] * coeff(tmp, x, L)
    end
    return p
end
function feynman_integral_branch_type(F::FeynmanIntegral, a::Vector{Int64},g::Vector{Int64}; aa=1,l=zeros(Int,nv(F.G)))
    ee = Edge.(F.G.edge)
    N = sum(a)
    f = signature_and_multiplicities(F.G, a)
    p = 0
    L = lis(F.G, N, l)
    x = F.S[2]  # Vector of elements x from the polynomial ring R
    q = F.S[3]  # Vector of elements q from the polynomial ring R
    z = F.S[4]  # Vector of elements z from the polynomial ring R
    sz=1
    g=2 .* g
    for k in eachindex(l)
        sz=sz*filter_term(inv_sfunction(z[k],aa),z,g)
    end
    sz=filter_term(sz,z,g)
    for i in eachindex(f)
        tmp=1
        for j in eachindex(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * filter_term(constterm(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])],aa,N),z,g)

                elseif f[i][2][j] == 0
                        tmp = tmp * filter_term(constterm(x[dst(ee[j])], x[src(ee[j])],z[src(ee[j])], z[dst(ee[j])],aa, N),z,g)

                elseif f[i][2][j] == -2
                        tmp = tmp *filter_term(loopterm(z[src(ee[j])],q[j],aa,a[j]),z,g)

                else 
                        tmp = tmp * filter_term(proterm(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])], q[j],f[i][2][j],aa, N),z,g)
                end 
        end
        p=p+f[i][1]*coeff(tmp,x,L)
       # p=filter_term(p,z,g)
        if i==length(f)
            pr=sz*p
            p=coeff(pr,z,g)
        end
    end
    return p
end
@doc raw"""
    feynman_integral_branch_type_order(  F::FeynmanGraph,   a::Vector{Int64}  ,Ω::Vector{Int64};l=zeros(Int,nv(F.G)))

 compute the Feynman Integral for a specified branch type `a` for a fixed ordering `Ω`
    
# Examples (without vertex contribution)

```julia
julia> G=FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> F=FeynmanIntegral(G);

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];

julia> feynman_integral_branch_type_order(F,a,Ω)
128*q[2]^2=4*q[3]^2*q[6]^2
```
     feynman_integral_branch_type_order( F::FeynmanGraph,   a::Vector{Int64}, Ω::Vector{Int64};aa=0,l=zeros(Int,nv(F.G)),g=zeros(Int,nv(F.G)))

# Examples (with vertex contribution)

```julia
julia> G=FeynmanGraph([(1, 2), (2, 3), (1, 3)])
FeynmanGraph([(1, 2), (2, 3), (1, 3)])

julia> F=FeynmanIntegral(G);

julia> a=[0,0,3];
julia> g=[1,0,0];
julia> Ω=[1,2,3];


julia> feynman_integral_branch_type_order(F,a,Ω,g)
115//6*q[3]^6
```
"""
function feynman_integral_branch_type_order(F::FeynmanIntegral, a::Vector{Int64}, o::Vector{Int64}; l=zeros(Int, nv(F.G)))
    ee = Edge.(F.G.edge)
    N = sum(a)
    f = signature_and_multiplicities_order(F.G, a,o)
    p = 0
    L = lis(F.G, N, l)
    
    x = F.S[2]  # Vector of elements x from the polynomial ring R
    q = F.S[3]  # Vector of elements q from the polynomial ring R


    for i in eachindex(f)
        tmp = 1  # Use the identity element of the polynomial ring R

        for j in eachindex(f[i][2])
            if f[i][2][j] == -1
                tmp *= constterm(x[src(ee[j])], x[dst(ee[j])], N)
            elseif f[i][2][j] == 0
                tmp *= constterm(x[dst(ee[j])], x[src(ee[j])], N)
            elseif f[i][2][j] == -2
                tmp *= loopterm(q[j], a[j])
            else
                tmp *= proterm(x[src(ee[j])], x[dst(ee[j])], q[j], f[i][2][j], N)
            end
        end

        p += f[i][1] * coeff(tmp, x, L)
    end
    return p
end
function feynman_integral_branch_type_order(F::FeynmanIntegral, a::Vector{Int64}, o::Vector{Int64}, g::Vector{Int64}; aa=1,l=zeros(Int,nv(F.G)))
    ee = Edge.(F.G.edge)
    N = sum(a)
    f = signature_and_multiplicities_order(F.G, a,o)
    p = 0
    L = lis(F.G, N, l)
    x = F.S[2]  # Vector of elements x from the polynomial ring R
    q = F.S[3]  # Vector of elements q from the polynomial ring R
    z = F.S[4]  # Vector of elements z from the polynomial ring R

    sz=1
    g=2 .* g

    for k in eachindex(l)
        sz=sz*filter_term(inv_sfunction(z[k],aa),z,g)
    end
    sz=filter_term(sz,z,g)
    for i in eachindex(f)
        tmp=1
        for j in eachindex(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * filter_term(constterm(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])],aa,N),z,g)

                elseif f[i][2][j] == 0
                        tmp = tmp * filter_term(constterm(x[dst(ee[j])], x[src(ee[j])],z[src(ee[j])], z[dst(ee[j])],aa, N),z,g)

                elseif f[i][2][j] == -2
                        tmp = tmp *filter_term(loopterm(z[src(ee[j])],q[j],aa,a[j]),z,g)

                else 
                        tmp = tmp * filter_term(proterm(x[src(ee[j])], x[dst(ee[j])], z[src(ee[j])], z[dst(ee[j])], q[j],f[i][2][j],aa, N),z,g)
                end 
        end
        p=p+f[i][1]*coeff(tmp,x,L)
       # p=filter_term(p,z,g)
        if i==length(f)
            pr=sz*p
            p=coeff(pr,z,g)
        end
    end
    return p
end
@doc raw"""
    feynman_integral_degree_order(  F::FeynmanGraph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(F.G)))

 compute the Feynman Integral for  all the partitions of the degree d  for a fixed ordering `Ω`
    
# Examples (without vertex contribution)
```julia
julia> G=FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> F=FeynmanIntegral(G);

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];

julia> feynman_integral_degree_order(F,Ω,4)
4*q[1]^4*q[2]^2*q[3]^2 + 4*q[1]^2*q[2]^4*q[5]^2 + 4*q[1]^2*q[2]^4*q[6]^2 + 4*q[1]^2*q[3]^4*q[5]^2 + 4*q[1]^2*q[3]^4*q[6]^2 + 176*q[2]^8 + 496*q[2]^6*q[3]^2 + 60*q[2]^6*q[5]^2 + 60*q[2]^6*q[6]^2 + 788*q[2]^4*q[3]^4 + 128*q[2]^4*q[3]^2*q[5]^2 + 128*q[2]^4*q[3]^2*q[6]^2 + 4*q[2]^4*q[4]^2*q[5]^2 + 4*q[2]^4*q[4]^2*q[6]^2 + 16*q[2]^4*q[5]^4 + 16*q[2]^4*q[6]^4 + 496*q[2]^2*q[3]^6 + 128*q[2]^2*q[3]^4*q[5]^2 + 128*q[2]^2*q[3]^4*q[6]^2 + 4*q[2]^2*q[3]^2*q[4]^4 + 48*q[2]^2*q[3]^2*q[5]^4 + 4*q[2]^2*q[3]^2*q[5]^2*q[6]^2 + 48*q[2]^2*q[3]^2*q[6]^4 + 176*q[3]^8 + 60*q[3]^6*q[5]^2 + 60*q[3]^6*q[6]^2 + 4*q[3]^4*q[4]^2*q[5]^2 + 4*q[3]^4*q[4]^2*q[6]^2 + 16*q[3]^4*q[5]^4 + 16*q[3]^4*q[6]^4
```
 
     feynman_integral_degree_order(  F::FeynmanGraph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(F.G)),g=zeros(Int,nv(F.G)))
 
# Examples (with vertex contribution)

```julia
julia> G=FeynmanGraph([(1, 2), (2, 3), (1, 3)])
 FeynmanGraph([(1, 2), (2, 3), (1, 3)])

julia> F=FeynmanIntegral(G);

julia> a=[0,0,3];
julia> g=[1,0,0];
julia> Ω=[1,2,3];

julia> feynman_integral_degree_order(F,Ω,3,g)
1//24*q[1]^4*q[2]^2 + 1//24*q[1]^4*q[3]^2 + 1//24*q[1]^2*q[2]^4 + 1//12*q[1]^2*q[2]^2*q[3]^2 + 1//24*q[1]^2*q[3]^4 + 1//24*q[2]^4*q[3]^2 + 1//24*q[2]^2*q[3]^4 + 115//6*q[3]^6
```
"""
function feynman_integral_degree_order( F::FeynmanIntegral,o::Vector{Int64},d::Integer,;l=zeros(Int,nv(F.G)) )
    # ve==F.G.edge
     indices = find_equal_pairs(F.G.edge)
     if isempty(indices)
         return feynman_integral_deg_order(F,o,d;l)
     else
         ee=Edge.(F.G.edge)
 
         re = Vector{Vector{Any}}()
         res = []
         L = partition(length(ee), d)
 
         while !isempty(L)
             ll = popfirst!(L)
 
             found = false
 
            
 
             if !found
                 ge = generate_permutation(ll, indices)
                 L = setdiff(L, ge)  # Remove the processed partition from L
                 if length(ge) == 1
                     kk= ge[1]
                     push!(res, feynman_integral_branch_type_order(F, kk,o;l))
                 else
                     k=ge[1]
                     f = feynman_integral_branch_type_order(F, k,o;l)
                     if f != 0
                         push!(res, f)
                         c1 = collect(coefficients(f))[1]
                         for i in eachindex(ge)[2:end]
                             li = ge[i]
                             r2 = c1 * vector_to_monomial(F,li)
                             push!(res,  r2)
 
                         end
                     end
                 end
             end
         end
         if isempty(res)
             return 0
         else
             return sum(res)
         end
     end
 end
function feynman_integral_degree_order( F::FeynmanIntegral,o::Vector{Int64},d::Integer,g::Vector{Int} ;aa=1,l=zeros(Int,nv(F.G)) )
    # ve==F.G.edge
     indices = find_equal_pairs(F.G.edge)
     if isempty(indices)
         return feynman_integral_deg_order(F,o,d,g;aa,l)
     else
         ee=Edge.(F.G.edge)
 
         re = Vector{Vector{Any}}()
         res = []
         L = partition(length(ee), d)
 
         while !isempty(L)
             ll = popfirst!(L)
 
             found = false
 
            
 
             if !found
                 ge = generate_permutation(ll, indices)
                 L = setdiff(L, ge)  # Remove the processed partition from L
                 if length(ge) == 1
                     kk= ge[1]
                     push!(res, feynman_integral_branch_type_order(F, kk,o,g;aa,l))
                 else
                     k=ge[1]
                     f = feynman_integral_branch_type_order(F, k,o,g;aa,l)
                     if f != 0
                         push!(res, f)
                         c1 = collect(coefficients(f))[1]
                         for i in eachindex(ge)[2:end]
                             li = ge[i]
                             r2 = c1 * vector_to_monomial(F,li)
                             push!(res,  r2)
 
                         end
                     end
                 end
             end
         end
         if isempty(res)
             return 0
         else
             return sum(res)
         end
     end
 end
function feynman_integral_deg_order(F::FeynmanIntegral,o::Vector{Int64},d::Integer;l=zeros(Int,nv(F.G)))
    ee=Edge.(F.G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type_order(F,ai,o;l) 
    end
    return sum
end 
function feynman_integral_deg_order(F::FeynmanIntegral,o::Vector{Int64},d::Integer,g;aa=1,l=zeros(Int,nv(F.G)))
    ee=Edge.(F.G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type_order(F,ai,o,g;aa,l) 
    end
    return sum
end

@doc raw"""
     feynman_integral_degree( F::FeynmanGraph,d::Integer;l=zeros(Int,nv(F.G)))

 compute the Feynman Integral for  over all the partitions of the degree d  for all  ordering `Ω`
    
# Examples (without vertex contribution)

```julia
julia> G=FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
 FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> feynman_integral_degree(F,3)

288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 + 32*q[1]^4*q[6]^2 + 8*q[1]^2*q[2]^2*q[5]^2 + 8*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[3]^2*q[5]^2 + 8*q[1]^2*q[3]^2*q[6]^2 + 24*q[2]^6 + 152*q[2]^4*q[3]^2 + 8*q[2]^4*q[5]^2 + 8*q[2]^4*q[6]^2 + 152*q[2]^2*q[3]^4 + 32*q[2]^2*q[3]^2*q[5]^2 + 32*q[2]^2*q[3]^2*q[6]^2 + 32*q[2]^2*q[4]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 8*q[2]^2*q[4]^2*q[6]^2 + 8*q[2]^2*q[5]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 8*q[2]^2*q[6]^4 + 24*q[3]^6 + 8*q[3]^4*q[5]^2 + 8*q[3]^4*q[6]^2 + 32*q[3]^2*q[4]^4 + 8*q[3]^2*q[4]^2*q[5]^2 + 8*q[3]^2*q[4]^2*q[6]^2 + 8*q[3]^2*q[5]^4 + 32*q[3]^2*q[5]^2*q[6]^2 + 8*q[3]^2*q[6]^4 + 288*q[4]^6 + 32*q[4]^4*q[5]^2 + 32*q[4]^4*q[6]^2 + 24*q[5]^6 + 152*q[5]^4*q[6]^2 + 152*q[5]^2*q[6]^4 + 24*q[6]^6
```
 --------------------------

     feynman_integral_degree(  F::FeynmanGraph,d::Integer;aa=0,l=zeros(Int,nv(F.G)),g=zeros(Int,nv(F.G)))
 
# Examples (with vertex contribution)
```julia
julia> G=FeynmanGraph([(1, 2), (2, 3), (1, 3)])
 FeynmanGraph([(1, 2), (2, 3), (1, 3)])

julia> F=FeynmanIntegral(G);

julia> a=[0,0,3];
julia> g=[1,0,0];

julia> feynman_integral_degree(F,3,g)
115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[3]^4 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 1//4*q[2]^2*q[3]^4 + 115//3*q[3]^6
```
 """
 function feynman_integral_degree(F::FeynmanIntegral, d::Int64 ; l=zeros(Int, nv(F.G)) )
    # ve==F.G.edge
     indices = find_equal_pairs(F.G.edge)
     if isempty(indices)
         return feynman_integral_deg(F,d;l)
     else
         ee=Edge.(F.G.edge)
 
         re = Vector{Vector{Any}}()
         res = []
         L = partition(length(ee), d)
 
         while !isempty(L)
             ll = popfirst!(L)
 
             found = false
             if !found
                 ge = generate_permutation(ll, indices)
                 L = setdiff(L, ge)  # Remove the processed partition from L
                 if length(ge) == 1
                     kk= ge[1]
                     push!(res, feynman_integral_branch_type(F, kk;l))
                 else
                     k=ge[1]
                     f = feynman_integral_branch_type(F, k;l)
                     if f != 0
                         push!(res, f)
                         c1 = collect(coefficients(f))[1]
                         for i in eachindex(ge)[2:end]
                             li = ge[i]
                             r2 = c1 * vector_to_monomial(F,li)
                             push!(res,  r2)
 
                         end
                     end
                 end
             end
         end
         if isempty(res)
             return 0
         else
             return sum(res)
         end
     end
 end
 function feynman_integral_degree( F::FeynmanIntegral,d::Integer,g::Vector{Int} ;aa=1,l=zeros(Int,nv(F.G)) )
    # ve==F.G.edge
     indices = find_equal_pairs(F.G.edge)
     if isempty(indices)
         return feynman_integral_deg(F,d,g;aa,l)
     else
         ee=Edge.(F.G.edge)
 
         re = Vector{Vector{Any}}()
         res = []
         L = partition(length(ee), d)
 
         while !isempty(L)
             ll = popfirst!(L)
 
             found = false
             if !found
                 ge = generate_permutation(ll, indices)
                 L = setdiff(L, ge)  # Remove the processed partition from L
                 if length(ge) == 1
                     kk= ge[1]
                     push!(res, feynman_integral_branch_type(F, kk,g;aa,l))
                 else
                     k=ge[1]
                     f = feynman_integral_branch_type(F, k,g;aa,l)
                     if f != 0
                         push!(res, f)
                         c1 = collect(coefficients(f))[1]
                         for i in eachindex(ge)[2:end]
                             li = ge[i]
                             r2 = c1 * vector_to_monomial(F,li)
                             push!(res,  r2)
 
                         end
                     end
                 end
             end
         end
         if isempty(res)
             return 0
         else
             return sum(res)
         end
     end
 end
 function feynman_integral_deg( F::FeynmanIntegral,d::Integer ;l=zeros(Int,nv(F.G)))
    ee=Edge.(F.G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type(F,ai;l) 
    end
    return sum
end 
function feynman_integral_deg( F::FeynmanIntegral,d::Integer,g ;aa=1,l=zeros(Int,nv(F.G)))
    ee=Edge.(F.G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type(F,ai,g;aa,l) 
    end
    return sum
end 

@doc raw"""
     feynman_integral_degree_sum_order( F::FeynmanGraph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(F.G)), g=zeros(Int, nv(F.G)))

compute the sum of all Feynman Integrals up to a certain degree d with a fixed ordering Ω

# Examples (without vertex contribution)

```julia
julia> G=FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
 FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];


julia> feynman_integral_degree_sum_order(F,Ω,3)
12*q[2]^6 + 76*q[2]^4*q[3]^2 + 4*q[2]^4*q[5]^2 + 4*q[2]^4*q[6]^2 + 76*q[2]^2*q[3]^4 + 16*q[2]^2*q[3]^2*q[5]^2 + 16*q[2]^2*q[3]^2*q[6]^2 + 4*q[2]^2*q[3]^2 + 12*q[3]^6 + 4*q[3]^4*q[5]^2 + 4*q[3]^4*q[6]^2
```
--------------------------

     feynman_integral_degree_sum_order( F::FeynmanGraph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(F.G)), g=zeros(Int, nv(F.G)))
 # Examples (with vertex contribution)

 ```julia
 julia> G=FeynmanGraph([(1, 2), (2, 3), (1, 3)])
 FeynmanGraph([(1, 2), (2, 3), (1, 3)])

julia> F=FeynmanIntegral(G);

julia> a=[0,0,3];
julia> g=[1,0,0];
julia> Ω=[1,2,3];


julia> feynman_integral_degree_sum_order(F,Ω,3,g)
1//24*q[1]^4*q[2]^2 + 1//24*q[1]^4*q[3]^2 + 1//24*q[1]^2*q[2]^4 + 1//12*q[1]^2*q[2]^2*q[3]^2 + 1//24*q[1]^2*q[2]^2 + 1//24*q[1]^2*q[3]^4 + 1//24*q[1]^2*q[3]^2 + 1//24*q[2]^4*q[3]^2 + 1//24*q[2]^2*q[3]^4 + 1//24*q[2]^2*q[3]^2 + 115//6*q[3]^6 + 19//8*q[3]^4 + 1//24*q[3]^2
```
 """
 function feynman_integral_degree_sum_order( F::FeynmanIntegral,o::Vector{Int64}, d::Union{Int64, Vector{Int64}};l=zeros(Int, nv(F.G)))
    res = 0
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree_order( F,o, i; l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree_order( F, o,i; l)
        end
    end
    return res
end
function feynman_integral_degree_sum_order( F::FeynmanIntegral,o::Vector{Int64}, d::Union{Int64, Vector{Int64}},g; aa=1,l=zeros(Int, nv(F.G)))
    res = 0
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree_order( F,o, i,g;aa, l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree_order( F, o,i,g;aa, l)
        end
    end
    return res
end
@doc raw"""
     feynman_integral_degree_sum( F::FeynmanGraph, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(F.G)))

 compute the sum of all Feynman Integrals up to a certain degree d for all  ordering Ω
    
# Examples (without vertex contribution)
```julia
julia> G=FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
FeynmanGraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> F=FeynmanIntegral(G);


julia> feynman_integral_degree_sum(F,3)
288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 + 32*q[1]^4*q[6]^2 + 8*q[1]^4 + 8*q[1]^2*q[2]^2*q[5]^2 + 8*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[3]^2*q[5]^2 + 8*q[1]^2*q[3]^2*q[6]^2 + 24*q[2]^6 + 152*q[2]^4*q[3]^2 + 8*q[2]^4*q[5]^2 + 8*q[2]^4*q[6]^2 + 152*q[2]^2*q[3]^4 + 32*q[2]^2*q[3]^2*q[5]^2 + 32*q[2]^2*q[3]^2*q[6]^2 + 8*q[2]^2*q[3]^2 + 32*q[2]^2*q[4]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 8*q[2]^2*q[4]^2*q[6]^2 + 8*q[2]^2*q[5]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 8*q[2]^2*q[6]^4 + 24*q[3]^6 + 8*q[3]^4*q[5]^2 + 8*q[3]^4*q[6]^2 + 32*q[3]^2*q[4]^4 + 8*q[3]^2*q[4]^2*q[5]^2 + 8*q[3]^2*q[4]^2*q[6]^2 + 8*q[3]^2*q[5]^4 + 32*q[3]^2*q[5]^2*q[6]^2 + 8*q[3]^2*q[6]^4 + 288*q[4]^6 + 32*q[4]^4*q[5]^2 + 32*q[4]^4*q[6]^2 + 8*q[4]^4 + 24*q[5]^6 + 152*q[5]^4*q[6]^2 + 152*q[5]^2*q[6]^4 + 8*q[5]^2*q[6]^2 + 24*q[6]^6
```
 --------------------------
 
     feynman_integral_degree_sum( F::FeynmanGraph, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(F.G)), g=zeros(Int, nv(F.G)))
 
# Examples (with vertex contribution)

```julia
julia> G=FeynmanGraph([(1, 2), (2, 3), (1, 3)])
 FeynmanGraph([(1, 2), (2, 3), (1, 3)])

julia> F=FeynmanIntegral(G);


julia> feynman_integral_degree_sum(F,3,g)
115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 19//4*q[1]^4 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[2]^2 + 1//4*q[1]^2*q[3]^4 + 1//4*q[1]^2*q[3]^2 + 1//12*q[1]^2 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 19//4*q[2]^4 + 1//4*q[2]^2*q[3]^4 + 1//4*q[2]^2*q[3]^2 + 1//12*q[2]^2 + 115//3*q[3]^6 + 19//4*q[3]^4 + 1//12*q[3]^2
```
 """
 function feynman_integral_degree_sum( F::FeynmanIntegral, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(F.G)))
    res = 0
    
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree( F, i; l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree( F, i; l)
        end
    end
    return res
end
function feynman_integral_degree_sum( F::FeynmanIntegral, d::Union{Int64, Vector{Int64}},g; aa=1,l=zeros(Int, nv(F.G)))
    res = 0
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree( F, i,g;aa, l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree( F, i,g;aa, l)
        end
    end
    return res
end
@doc raw"""
     sum_of_coeff(p::QQMPolyRingElem)

compute the sum of coefficient of the polynomial p. 

```julia
julia> f=3*x[1]^6 + 2*x[1]^5*x[2] + x[1]^4*x[2]^2

julia> sum_of_coeff(f)
 6
```
 """
function sum_of_coeff(p::QQMPolyRingElem)
    coeffs_dict = coefficients(p)
    coeffs_array = collect(values(coeffs_dict))
    return sum(coeffs_array)
end
@doc raw"""
     substitute(p::Union{QQMPolyRingElem, Int64})

 replace all the variables by the first variable of p. 
 With `x=\[x_1,x2,x3 \]` and `p(x_1,x2,x3)`, `substitute(x,p)` returns `p(x1,x_1,x_1)`

```julia
julia> f=x[1]*x[2]+x[1]^3*x[2]+5x[1]^6-2x[3]*x[2]
  5*x[1]^6 + x[1]^3*x[2] + x[1]*x[2] - 2*x[2]*x[3]

julia> substitute(f)
 5*x[1]^6 + x[1]^4 - x[1]^2
```
 """
 function substitute(p::Union{QQMPolyRingElem, Int64})
    if typeof(p)==Int64 || p==zero(p)
        return zero(p)
    else
    q=vars(p)
        T=parent(p)
        r=gens(T)       
        s=findfirst(isequal(q[1]),r)
        r[s:end] .= q[1]
    end
    return evaluate(p,r)
end
