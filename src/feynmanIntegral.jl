@doc raw"""
    feynman_integral_branch_type( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;l=zeros(Int,nv(G)))

 compute the Feynman Integral for a specified branch type `a` for all ordering `Ω`
    
# Examples (without vertex contribution)
 
```julia
julia> G=graph(ve)
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> feynman_integral_branch_type(x,q,G,a)
256*q[2]^2*q[3]*q[6]
```

     feynman_integral_branch_type(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))

# Examples (with vertex contribution)

```julia
julia> G=graph(ve)
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> feynman_integral_branch_type(x,q,z,G,a,aa=1,g=[1,0,0])
 115//3*q[3]^3
```
"""
function feynman_integral_branch_type( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities(G, a)
    p=0
    L=lis(G,N,l)
    for i in 1:length(f)
        tmp=1
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
function feynman_integral_branch_type(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities(G, a)
    L=lis(G,N,l)
    g=2 .* g
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
julia> G=graph(ve)
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];

julia> feynman_integral_branch_type_order(x,q,G,a,Ω)
128*q[2]^2*q[3]*q[6]
```
     feynman_integral_branch_type_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}, Ω::Vector{Int64};aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))

# Examples (with vertex contribution)

```julia
julia> G=graph(ve)
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> Ω=[1,2,3];

julia> feynman_integral_branch_type_order(x,q,z,G,a,Ω,aa=1,g=[1,0,0])
 115//6*q[3]^3
```
"""
function feynman_integral_branch_type_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}  ,o::Vector{Int64};l=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities_order(G, a,o)
    p=parent(x[1])(0) 
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
function feynman_integral_branch_type_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64},o::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee = Edge.(G.edge)
    N = sum(a)
    f = signature_and_multiplicities_order(G, a,o)
    L=lis(G,N,l)
    g=2 .* g
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
    feynman_integral_degree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)))

 compute the Feynman Integral for  all the partitions of the degree d  for a fixed ordering `Ω`
    
# Examples (without vertex contribution)
```julia
julia> G=graph(ve)

 graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];

julia> feynman_integral_degree_order(x,q,G,Ω,4)

 4*q[1]^2*q[2]*q[3] + 4*q[1]*q[2]^2*q[5] + 4*q[1]*q[2]^2*q[6] + 4*q[1]*q[3]^2*q[5] + 4*q[1]*q[3]^2*q[6] + 176*q[2]^4 + 496*q[2]^3*q[3] + 60*q[2]^3*q[5] + 60*q[2]^3*q[6] + 788*q[2]^2*q[3]^2 + 128*q[2]^2*q[3]*q[5] + 128*q[2]^2*q[3]*q[6] + 4*q[2]^2*q[4]*q[5] + 4*q[2]^2*q[4]*q[6] + 16*q[2]^2*q[5]^2 + 16*q[2]^2*q[6]^2 + 496*q[2]*q[3]^3 + 128*q[2]*q[3]^2*q[5] + 128*q[2]*q[3]^2*q[6] + 4*q[2]*q[3]*q[4]^2 + 48*q[2]*q[3]*q[5]^2 + 4*q[2]*q[3]*q[5]*q[6] + 48*q[2]*q[3]*q[6]^2 + 176*q[3]^4 + 60*q[3]^3*q[5] + 60*q[3]^3*q[6] + 4*q[3]^2*q[4]*q[5] + 4*q[3]^2*q[4]*q[6] + 16*q[3]^2*q[5]^2 + 16*q[3]^2*q[6]^2
```
 
     feynman_integral_degree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
 
# Examples (with vertex contribution)

```julia
julia> G=graph(ve)
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> Ω=[1,2,3];

julia> feynman_integral_degree_order(x,q,z,G,o,3,aa=1,g=[1,0,0])
 1//24*q[1]^2*q[2] + 1//24*q[1]^2*q[3] + 1//24*q[1]*q[2]^2 + 1//12*q[1]*q[2]*q[3] + 1//24*q[1]*q[3]^2 + 1//24*q[2]^2*q[3] + 1//24*q[2]*q[3]^2 + 115//6*q[3]^3
```
"""
function feynman_integral_degree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type_order(x,q,G,ai,o;l) 
    end
    return sum
end 
function feynman_integral_degree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type_order(x,q,z,G,ai,o;aa,l,g) 
    end
    return sum
end 

@doc raw"""
     feynman_integral_degree(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,d::Integer;l=zeros(Int,nv(G)))

 compute the Feynman Integral for  over all the partitions of the degree d  for all  ordering `Ω`
    
# Examples (without vertex contribution)

```julia
julia> G=graph(ve)

 graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> feynman_integral_degree(x,q,G,3)

 288*q[1]^3 + 32*q[1]^2*q[2] + 32*q[1]^2*q[3] + 32*q[1]^2*q[5] + 32*q[1]^2*q[6] + 8*q[1]*q[2]*q[5] + 8*q[1]*q[2]*q[6] + 8*q[1]*q[3]*q[5] + 8*q[1]*q[3]*q[6] + 24*q[2]^3 + 152*q[2]^2*q[3] + 8*q[2]^2*q[5] + 8*q[2]^2*q[6] + 152*q[2]*q[3]^2 + 32*q[2]*q[3]*q[5] + 32*q[2]*q[3]*q[6] + 32*q[2]*q[4]^2 + 8*q[2]*q[4]*q[5] + 8*q[2]*q[4]*q[6] + 8*q[2]*q[5]^2 + 32*q[2]*q[5]*q[6] + 8*q[2]*q[6]^2 + 24*q[3]^3 + 8*q[3]^2*q[5] + 8*q[3]^2*q[6] + 32*q[3]*q[4]^2 + 8*q[3]*q[4]*q[5] + 8*q[3]*q[4]*q[6] + 8*q[3]*q[5]^2 + 32*q[3]*q[5]*q[6] + 8*q[3]*q[6]^2 + 288*q[4]^3 + 32*q[4]^2*q[5] + 32*q[4]^2*q[6] + 24*q[5]^3 + 152*q[5]^2*q[6] + 152*q[5]*q[6]^2 + 24*q[6]^3
```
 --------------------------

     feynman_integral_degree( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
 
# Examples (with vertex contribution)
```julia
julia> G=graph(ve)
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> feynman_integral_degree(x,q,z,G,3,aa=1,g=[1,0,0])
 115//3*q[1]^3 + 1//4*q[1]^2*q[2] + 1//4*q[1]^2*q[3] + 1//4*q[1]*q[2]^2 + 1//2*q[1]*q[2]*q[3] + 1//4*q[1]*q[3]^2 + 115//3*q[2]^3 + 1//4*q[2]^2*q[3] + 1//4*q[2]*q[3]^2 + 115//3*q[3]^3
```
 """
function feynman_integral_degree(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,d::Integer;l=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type(x,q,G,ai;l) 
    end
    return sum
end 
function feynman_integral_degree( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))
    ee=Edge.(G.edge)
    a=partition(length(ee),d) 
    sum=0
    for ai in a
         sum=sum+feynman_integral_branch_type(x,q,z,G,ai;aa,l,g) 
    end
    return sum
end 

@doc raw"""
     feynman_integral_degree_sum_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))

compute the sum of all Feynman Integrals up to a certain degree d with a fixed ordering Ω

# Examples (without vertex contribution)

```julia
julia> G=graph(ve)
 graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

julia> a=[0,2,1,0,0,1];

julia> Ω=[1,3,4,2];

julia> feynman_integral_degree_sum_order(x,q,G,Ω,3)
 12*q[2]^3 + 76*q[2]^2*q[3] + 4*q[2]^2*q[5] + 4*q[2]^2*q[6] + 76*q[2]*q[3]^2 + 16*q[2]*q[3]*q[5] + 16*q[2]*q[3]*q[6] + 4*q[2]*q[3] + 12*q[3]^3 + 4*q[3]^2*q[5] + 4*q[3]^2*q[6]
```
--------------------------

     feynman_integral_degree_sum_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))
 # Examples (with vertex contribution)

 ```julia
 julia> G=graph(ve)
 graph([(1, 2), (2, 3), (1, 3)])

julia> a=[0,0,3];

julia> Ω=[1,2,3];

julia> feynman_integral_degree_sum_order(x,q,z,G,o,3,aa=1,g=[1,0,0])

 1//24*q[1]^2*q[2] + 1//24*q[1]^2*q[3] + 1//24*q[1]*q[2]^2 + 1//12*q[1]*q[2]*q[3] + 1//24*q[1]*q[2] + 1//24*q[1]*q[3]^2 + 1//24*q[1]*q[3] + 1//24*q[2]^2*q[3] + 1//24*q[2]*q[3]^2 + 1//24*q[2]*q[3] + 115//6*q[3]^3 + 19//8*q[3]^2 + 1//24*q[3]
```
 """
function feynman_integral_degree_sum_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))
    res = zero(x[1])
    
    if typeof(d) <: Integer
        for i in 1:d
            res+=feynman_integral_degree_order(x,q,G,o,i;l)
        end
    else
        for i in d[1]:d[end]
            res+=feynman_integral_degree_order(x,q,G,o,i;l)
        end
    end
    
    return res
end
function feynman_integral_degree_sum_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))
    res = zero(x[1])
    
    if typeof(d) <: Integer
        for i in 1:d
            res+=feynman_integral_degree_order(x,q,z,G,o,i;aa,l,g)
        end
    else
        for i in d[1]:d[end]
            res+=feynman_integral_degree_order(x,q,z,G,o,i;aa,l,g)
        end
    end
    
    return res
end
@doc raw"""
     feynman_integral_degree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(G)))

 compute the sum of all Feynman Integrals up to a certain degree d for all  ordering Ω
    
# Examples (without vertex contribution)
```julia
julia> G=graph(ve)
graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])

 julia> feynman_integral_degree_sum(x,q,G,3)

 288*q[1]^3 + 32*q[1]^2*q[2] + 32*q[1]^2*q[3] + 32*q[1]^2*q[5] + 32*q[1]^2*q[6] + 8*q[1]^2 + 8*q[1]*q[2]*q[5] + 8*q[1]*q[2]*q[6] + 8*q[1]*q[3]*q[5] + 8*q[1]*q[3]*q[6] + 24*q[2]^3 + 152*q[2]^2*q[3] + 8*q[2]^2*q[5] + 8*q[2]^2*q[6] + 152*q[2]*q[3]^2 + 32*q[2]*q[3]*q[5] + 32*q[2]*q[3]*q[6] + 8*q[2]*q[3] + 32*q[2]*q[4]^2 + 8*q[2]*q[4]*q[5] + 8*q[2]*q[4]*q[6] + 8*q[2]*q[5]^2 + 32*q[2]*q[5]*q[6] + 8*q[2]*q[6]^2 + 24*q[3]^3 + 8*q[3]^2*q[5] + 8*q[3]^2*q[6] + 32*q[3]*q[4]^2 + 8*q[3]*q[4]*q[5] + 8*q[3]*q[4]*q[6] + 8*q[3]*q[5]^2 + 32*q[3]*q[5]*q[6] + 8*q[3]*q[6]^2 + 288*q[4]^3 + 32*q[4]^2*q[5] + 32*q[4]^2*q[6] + 8*q[4]^2 + 24*q[5]^3 + 152*q[5]^2*q[6] + 152*q[5]*q[6]^2 + 8*q[5]*q[6] + 24*q[6]^3
```
 --------------------------
 
     feynman_integral_degree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))
 
# Examples (with vertex contribution)

```julia
julia> G=graph(ve)
 graph([(1, 2), (2, 3), (1, 3)])

julia> feynman_integral_degree_sum(x,q,z,G,3,aa=1,g=[1,0,0])

 115//3*q[1]^3 + 1//4*q[1]^2*q[2] + 1//4*q[1]^2*q[3] + 19//4*q[1]^2 + 1//4*q[1]*q[2]^2 + 1//2*q[1]*q[2]*q[3] + 1//4*q[1]*q[2] + 1//4*q[1]*q[3]^2 + 1//4*q[1]*q[3] + 1//12*q[1] + 115//3*q[2]^3 + 1//4*q[2]^2*q[3] + 19//4*q[2]^2 + 1//4*q[2]*q[3]^2 + 1//4*q[2]*q[3] + 1//12*q[2] + 115//3*q[3]^3 + 19//4*q[3]^2 + 1//12*q[3]
```
 """
function feynman_integral_degree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(G)))
    res = zero(x[1])
    
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree(x, q, G, i; l)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree(x, q, G, i; l)
        end
    end
    
    return res
end

function feynman_integral_degree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))
    res = zero(x[1])
    
    if typeof(d) <: Integer
        for i in 1:d
            res += feynman_integral_degree(x, q,z, G, i; aa, l, g)
        end
    else
        for i in d[1]:d[end]
            res += feynman_integral_degree(x, q ,z, G, i; aa, l, g)
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
julia> f=x[1]*x[2]+x[1]^3*x[2]+5x[1]^6-2x[3]*x[2]
  5*x[1]^6 + x[1]^3*x[2] + x[1]*x[2] - 2*x[2]*x[3]

julia> substitute(x,f)
 5*x[1]^6 + x[1]^4 - x[1]^2
```
 """
function substitute(q::Vector{QQMPolyRingElem},p::Union{QQMPolyRingElem, Int64})
    T=parent(q[1])
    if typeof(p)==Int64
        return 0
    else
        r=gens(T)
        s=findfirst(isequal(q[1]),r)
        r[s:end] .= q[1]
    end
    return evaluate(p,r)
end
