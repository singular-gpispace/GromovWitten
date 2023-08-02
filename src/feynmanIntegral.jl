@doc raw"""
    feynman_integral_branch_type( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;l=zeros(Int,nv(G)))

 compute the Feynman Integral for a specified branch type `a` for all ordering `Ω`
    
# Examples (without vertex contribution)
 
```julia
julia> G=graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];
julia> feynman_integral_branch_type(G,a)
256*q[2]^4*q[3]^2*q[6]^2
```

     feynman_integral_branch_type(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))

# Examples (with vertex contribution)

```julia
julia> G=graph([(1, 2), (2, 3), (1, 3)])
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];
julia> g=[1,0,0];
julia> feynman_integral_branch_type(G,a,g)
 115//3*q[3]^6
```
"""
function feynman_integral_branch_type(G::graph, a::Vector{Int64};l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities(G, a)
    p = 0
    L = lis(G, N,l)
    S=@polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)],z[1:nv(G)])
    x = S[2]  # Vector of elements x from the polynomial ring R
    q = S[3]  # Vector of elements q from the polynomial ring R

    for i in 1:length(f)
        tmp = 1  # Use the identity element of the polynomial ring R

        for j in 1:length(f[i][2])
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

function feynman_integral_branch_type(  G::graph ,a::Vector{Int64}, g::Vector{Int64}; aa=1,l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities(G, a)
    L=lis(G,N,l)
    g=2 .* g
    
    S=@polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)],z[1:nv(G)])
    x = S[2]  # Vector of elements x from the polynomial ring R
    q = S[3]  # Vector of elements q from the polynomial ring R
    z = S[4]  # Vector of elements q from the polynomial ring R
    p=0
     sz=1
    for k in 1:length(l)
        sz=sz*filter_term(inv_sfunction(z[k],aa),z,g)
    end
    sz=filter_term(sz,z,g)
    for i in 1:length(f)
        tmp=1
        for j in 1:length(f[i][2])

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
    feynman_integral_branch_type_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}  ,Ω::Vector{Int64};l=zeros(Int,nv(G)))

 compute the Feynman Integral for a specified branch type `a` for a fixed ordering `Ω`
    
# Examples (without vertex contribution)

```julia
julia> G=graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];
julia> feynman_integral_branch_type_order(G,a,Ω)
128*q[2]^2=4*q[3]^2*q[6]^2
```
     feynman_integral_branch_type_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}, Ω::Vector{Int64};aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))

# Examples (with vertex contribution)

```julia
julia> G=graph([(1, 2), (2, 3), (1, 3)])
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> Ω=[1,2,3];

julia> feynman_integral_branch_type_order(G,a,Ω,g)
 115//6*q[3]^6
```
"""
function feynman_integral_branch_type_order(  G::graph ,a::Vector{Int64}  ,o::Vector{Int64};l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities_order(G, a,o)
    S=@polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)] )
    x = S[2]  # Vector of elements x from the polynomial ring R
    q = S[3]  # Vector of elements q from the polynomial ring R
    p=S[1](0) 
    L=lis(G,N,l)
    
    for i in 1:length(f)
        tmp=1
        tm=1
        for j in 1:length(f[i][2])

                if f[i][2][j] == -1
                        tmp = tmp * constterm(x[src(ee[j])], x[dst(ee[j])],N)

                elseif f[i][2][j] == 0
                        tmp = tmp * constterm(x[dst(ee[j])], x[src(ee[j])], N)
                elseif f[i][2][j] == -2
                        tmp = tmp *loopterm(q[j],a[j])

                else 
                
                        tmp = tmp * proterm(x[src(ee[j])], x[dst(ee[j])], q[j],f[i][2][j], N)

                end 
        end
            #p=p+f[i][1]*tm
       # tmp=filter_term(tmp,x,L)

        p=p+f[i][1]*coeff(tmp,x,L)
       
    end
    return p
end
function feynman_integral_branch_type_order( G::graph ,a::Vector{Int64},o::Vector{Int64},g;aa=1,l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities_order(G, a,o)
    L=lis(G,N,l)
    S=@polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)],z[1:nv(G)])
    x = S[2]  # Vector of elements x from the polynomial ring R
    q = S[3]  # Vector of elements q from the polynomial ring R
    z = S[4]  # Vector of elements z from the polynomial ring R
    p=S[1](0) 
    g=2 .* g
     sz=1
    for k in 1:length(l)
        sz=sz*filter_term(inv_sfunction(z[k],aa),z,g)
    end
    sz=filter_term(sz,z,g)
    for i in 1:length(f)
        tmp=1
        for j in 1:length(f[i][2])

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
    feynman_integral_degree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)))

 compute the Feynman Integral for  all the partitions of the degree d  for a fixed ordering `Ω`
    
# Examples (without vertex contribution)
```julia
julia> G=graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

 graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];
julia> feynman_integral_degree_order(G,Ω,4)
4*q[1]^4*q[2]^2*q[3]^2 + 4*q[1]^2*q[2]^4*q[5]^2 + 4*q[1]^2*q[2]^4*q[6]^2 + 4*q[1]^2*q[3]^4*q[5]^2 + 4*q[1]^2*q[3]^4*q[6]^2 + 176*q[2]^8 + 496*q[2]^6*q[3]^2 + 60*q[2]^6*q[5]^2 + 60*q[2]^6*q[6]^2 + 788*q[2]^4*q[3]^4 + 128*q[2]^4*q[3]^2*q[5]^2 + 128*q[2]^4*q[3]^2*q[6]^2 + 4*q[2]^4*q[4]^2*q[5]^2 + 4*q[2]^4*q[4]^2*q[6]^2 + 16*q[2]^4*q[5]^4 + 16*q[2]^4*q[6]^4 + 496*q[2]^2*q[3]^6 + 128*q[2]^2*q[3]^4*q[5]^2 + 128*q[2]^2*q[3]^4*q[6]^2 + 4*q[2]^2*q[3]^2*q[4]^4 + 48*q[2]^2*q[3]^2*q[5]^4 + 4*q[2]^2*q[3]^2*q[5]^2*q[6]^2 + 48*q[2]^2*q[3]^2*q[6]^4 + 176*q[3]^8 + 60*q[3]^6*q[5]^2 + 60*q[3]^6*q[6]^2 + 4*q[3]^4*q[4]^2*q[5]^2 + 4*q[3]^4*q[4]^2*q[6]^2 + 16*q[3]^4*q[5]^4 + 16*q[3]^4*q[6]^4
```
 
     feynman_integral_degree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
 
# Examples (with vertex contribution)

```julia
julia> G=graph([(1, 2), (2, 3), (1, 3)])
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> Ω=[1,2,3];
julia> feynman_integral_degree_order(G,Ω,3,g)
1//24*q[1]^4*q[2]^2 + 1//24*q[1]^4*q[3]^2 + 1//24*q[1]^2*q[2]^4 + 1//12*q[1]^2*q[2]^2*q[3]^2 + 1//24*q[1]^2*q[3]^4 + 1//24*q[2]^4*q[3]^2 + 1//24*q[2]^2*q[3]^4 + 115//6*q[3]^6
```
"""
function feynman_integral_degree_order(G::graph,o::Vector{Int64},d::Integer;l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type_order(G,ai,o;l) 
    end
    return sum
end 
function feynman_integral_degree_order(G::graph,o::Vector{Int64},d::Integer,g;aa=1,l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type_order(G,ai,o,g;aa,l) 
    end
    return sum
end 

@doc raw"""
     feynman_integral_degree(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,d::Integer;l=zeros(Int,nv(G)))

 compute the Feynman Integral for  over all the partitions of the degree d  for all  ordering `Ω`
    
# Examples (without vertex contribution)

```julia
julia> G=graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
 graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];
julia> feynman_integral_degree(G,3)
288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 + 32*q[1]^4*q[6]^2 + 8*q[1]^2*q[2]^2*q[5]^2 + 8*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[3]^2*q[5]^2 + 8*q[1]^2*q[3]^2*q[6]^2 + 24*q[2]^6 + 152*q[2]^4*q[3]^2 + 8*q[2]^4*q[5]^2 + 8*q[2]^4*q[6]^2 + 152*q[2]^2*q[3]^4 + 32*q[2]^2*q[3]^2*q[5]^2 + 32*q[2]^2*q[3]^2*q[6]^2 + 32*q[2]^2*q[4]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 8*q[2]^2*q[4]^2*q[6]^2 + 8*q[2]^2*q[5]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 8*q[2]^2*q[6]^4 + 24*q[3]^6 + 8*q[3]^4*q[5]^2 + 8*q[3]^4*q[6]^2 + 32*q[3]^2*q[4]^4 + 8*q[3]^2*q[4]^2*q[5]^2 + 8*q[3]^2*q[4]^2*q[6]^2 + 8*q[3]^2*q[5]^4 + 32*q[3]^2*q[5]^2*q[6]^2 + 8*q[3]^2*q[6]^4 + 288*q[4]^6 + 32*q[4]^4*q[5]^2 + 32*q[4]^4*q[6]^2 + 24*q[5]^6 + 152*q[5]^4*q[6]^2 + 152*q[5]^2*q[6]^4 + 24*q[6]^6
```
 --------------------------

     feynman_integral_degree( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
 
# Examples (with vertex contribution)
```julia
julia> G=graph([(1, 2), (2, 3), (1, 3)])
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];
julia> feynman_integral_degree(G,3,g)
115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[3]^4 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 1//4*q[2]^2*q[3]^4 + 115//3*q[3]^6
```
 """
 function feynman_integral_degree( G::graph,d::Integer ;l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type(G,ai;l) 
    end
    return sum
end 
function feynman_integral_degree( G::graph,d::Integer,g ;aa=1,l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type(G,ai,g;aa,l) 
    end
    return sum
end 

@doc raw"""
     feynman_integral_degree_sum_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))

compute the sum of all Feynman Integrals up to a certain degree d with a fixed ordering Ω

# Examples (without vertex contribution)

```julia
julia> G=graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
 graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];

julia> feynman_integral_degree_sum_order(G,3,Ω)
12*q[2]^6 + 76*q[2]^4*q[3]^2 + 4*q[2]^4*q[5]^2 + 4*q[2]^4*q[6]^2 + 76*q[2]^2*q[3]^4 + 16*q[2]^2*q[3]^2*q[5]^2 + 16*q[2]^2*q[3]^2*q[6]^2 + 4*q[2]^2*q[3]^2 + 12*q[3]^6 + 4*q[3]^4*q[5]^2 + 4*q[3]^4*q[6]^2
```
--------------------------

     feynman_integral_degree_sum_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))
 # Examples (with vertex contribution)

 ```julia
 julia> G=graph([(1, 2), (2, 3), (1, 3)])
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> Ω=[1,2,3];

julia> feynman_integral_degree_sum_order(G,3,Ω,g)
1//24*q[1]^4*q[2]^2 + 1//24*q[1]^4*q[3]^2 + 1//24*q[1]^2*q[2]^4 + 1//12*q[1]^2*q[2]^2*q[3]^2 + 1//24*q[1]^2*q[2]^2 + 1//24*q[1]^2*q[3]^4 + 1//24*q[1]^2*q[3]^2 + 1//24*q[2]^4*q[3]^2 + 1//24*q[2]^2*q[3]^4 + 1//24*q[2]^2*q[3]^2 + 115//6*q[3]^6 + 19//8*q[3]^4 + 1//24*q[3]^2
```
 """
 function feynman_integral_degree_sum_order( G::graph, d::Union{Int64, Vector{Int64}},o::Vector{Int64};l=zeros(Int, nv(G)))
    res = 0
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree_order( G,o, i; l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree_order( G, o,i; l)
        end
    end
    return res
end
function feynman_integral_degree_sum_order( G::graph, d::Union{Int64, Vector{Int64}},o::Vector{Int64},g; aa=1,l=zeros(Int, nv(G)))
    res = 0
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree_order( G,o, i,g;aa, l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree_order( G, o,i,g;aa, l)
        end
    end
    return res
end
@doc raw"""
     feynman_integral_degree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(G)))

 compute the sum of all Feynman Integrals up to a certain degree d for all  ordering Ω
    
# Examples (without vertex contribution)
```julia
julia> G=graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> feynman_integral_degree_sum(G,3)
288*q[1]^6 + 32*q[1]^4*q[2]^2 + 32*q[1]^4*q[3]^2 + 32*q[1]^4*q[5]^2 + 32*q[1]^4*q[6]^2 + 8*q[1]^4 + 8*q[1]^2*q[2]^2*q[5]^2 + 8*q[1]^2*q[2]^2*q[6]^2 + 8*q[1]^2*q[3]^2*q[5]^2 + 8*q[1]^2*q[3]^2*q[6]^2 + 24*q[2]^6 + 152*q[2]^4*q[3]^2 + 8*q[2]^4*q[5]^2 + 8*q[2]^4*q[6]^2 + 152*q[2]^2*q[3]^4 + 32*q[2]^2*q[3]^2*q[5]^2 + 32*q[2]^2*q[3]^2*q[6]^2 + 8*q[2]^2*q[3]^2 + 32*q[2]^2*q[4]^4 + 8*q[2]^2*q[4]^2*q[5]^2 + 8*q[2]^2*q[4]^2*q[6]^2 + 8*q[2]^2*q[5]^4 + 32*q[2]^2*q[5]^2*q[6]^2 + 8*q[2]^2*q[6]^4 + 24*q[3]^6 + 8*q[3]^4*q[5]^2 + 8*q[3]^4*q[6]^2 + 32*q[3]^2*q[4]^4 + 8*q[3]^2*q[4]^2*q[5]^2 + 8*q[3]^2*q[4]^2*q[6]^2 + 8*q[3]^2*q[5]^4 + 32*q[3]^2*q[5]^2*q[6]^2 + 8*q[3]^2*q[6]^4 + 288*q[4]^6 + 32*q[4]^4*q[5]^2 + 32*q[4]^4*q[6]^2 + 8*q[4]^4 + 24*q[5]^6 + 152*q[5]^4*q[6]^2 + 152*q[5]^2*q[6]^4 + 8*q[5]^2*q[6]^2 + 24*q[6]^6
```
 --------------------------
 
     feynman_integral_degree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))
 
# Examples (with vertex contribution)

```julia
julia> G=graph([(1, 2), (2, 3), (1, 3)])
 graph([(1, 2), (2, 3), (1, 3)])

julia> feynman_integral_degree_sum(G,3,g)
115//3*q[1]^6 + 1//4*q[1]^4*q[2]^2 + 1//4*q[1]^4*q[3]^2 + 19//4*q[1]^4 + 1//4*q[1]^2*q[2]^4 + 1//2*q[1]^2*q[2]^2*q[3]^2 + 1//4*q[1]^2*q[2]^2 + 1//4*q[1]^2*q[3]^4 + 1//4*q[1]^2*q[3]^2 + 1//12*q[1]^2 + 115//3*q[2]^6 + 1//4*q[2]^4*q[3]^2 + 19//4*q[2]^4 + 1//4*q[2]^2*q[3]^4 + 1//4*q[2]^2*q[3]^2 + 1//12*q[2]^2 + 115//3*q[3]^6 + 19//4*q[3]^4 + 1//12*q[3]^2
```
 """
 function feynman_integral_degree_sum( G::graph, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(G)))
    res = 0
    
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree( G, i; l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree( G, i; l)
        end
    end
    return res
end
function feynman_integral_degree_sum( G::graph, d::Union{Int64, Vector{Int64}},g; aa=1,l=zeros(Int, nv(G)))
    res = 0
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree( G, i,g;aa, l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree( G, i,g;aa, l)
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
     substitute(q::Vector{QQMPolyRingElem},p::Union{QQMPolyRingElem, Int64})

 replace all the variables by the first variable of p. 
 With `x=\[x_1,x2,x3 \]` and `p(x_1,x2,x3)`, `substitute(x,p)` returns `p(x1,x_1,x_1)`

```julia
julia> R,x=@polynomial_ring(QQ,x[1:4])

julia> f=x[1]*x[2]+x[1]^3*x[2]+5x[1]^6-2x[3]*x[2]
  5*x[1]^6 + x[1]^3*x[2] + x[1]*x[2] - 2*x[2]*x[3]

julia> substitute(f)
 5*x[1]^6 + x[1]^4 - x[1]^2
```
 """
 function substitute(p::Union{fmpq_mpoly, Int64})
    q=vars(p)
    if typeof(p)==Int64
        return 0
    else
        T=parent(p)
        r=gens(T)       
        s=findfirst(isequal(q[1]),r)
        r[s:end] .= q[1]
    end
    return evaluate(p,r)
end 