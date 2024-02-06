###############################################################################
#                                                                             #
#     coefterm.jl : compute differents  coefficients of polynomial.           #
#                                                                             #
###############################################################################

#  coefterm compute the coefficients term x1^0,...,xn^0, By I(q)=coeff[x1^0,...,xn^0](P(x,q))
#Where I(q) is the Feynman Integral and P(x,q) the Propagator.
#=function coefterm(x::Vector,q::Vector,G::FeynmanGraph ,p::QQMPolyRingElem,d::Integer;l=zeros(Int,nv(G)))
   #here l is leak vector of the graph G.
    ee=edges(G) 
    G=DiGraph(edges(G)) # convert from graph to Graph (so we can use nv(G))
    L=zeros(Int,nv(G)) #  define zeros Vector of length nv(G)
    for ev in ee
        if src(ev) == dst(ev) #checking for loop term. 

        else
            L[src(ev)]=L[src(ev)]+d
            L[dst(ev)]=L[dst(ev)]+d
        end
    end # So L =d*ei where ei is the number of valence of each vertex
    L=L .+l # adding l to the vector L. 
    p=coeff(p,x,L) # compute the coefficients of degree x1^l1,...,xn^ln.
    return p
end =#

function replace(vector)
    result_vector = [x == -1 ? 0 : (x == 0 ? -1 : x) for x in vector]
    return result_vector
end
function count_zero(arr)
    count_zeros = 0

    for element in arr
        if element == 0
            count_zeros += 1
        end
    end
    return count_zeros
end
@doc raw"""
   partition(k::Integer, n::Integer) 


#Examples

This function returns the number of partitions of $n$ into fixed  $k$ parts. 



```jldoctest Gromov
julia> using GromovWitten

julia> partition(3, 4)
15-element Vector{Vector{Int64}}:
 [4, 0, 0]
 [3, 1, 0]
 [3, 0, 1]
 [2, 2, 0]
 [2, 1, 1]
 [2, 0, 2]
 [1, 3, 0]
 [1, 2, 1]
 [1, 1, 2]
 [1, 0, 3]
 [0, 4, 0]
 [0, 3, 1]
 [0, 2, 2]
 [0, 1, 3]
 [0, 0, 4]
```
"""
function partition(n::Integer,k::Integer)
    if(k==0)
        return [[0]]
    end
    if(n==0)
        return [[0] *k]
    end
    p=with_replacement_combinations(1:n,k)
    return map(A -> [sum(A .== i) for i in 1:n],p)
end
   
#give the position of the vertices xi in the list L. 
function preimg(L::Vector{Int64}, xi::Int64)
    for (i,Li) in enumerate(L)
        if Li==xi
            return i # i is then the position of xi in L
        end
    end
end

@doc raw"""
    flip_signature(F::FeynmanGraph ,p::Vector{Int64},a::Vector{Int64})

Let   Ω=[x1,...,xn] be  a given Order and $a$  a branche type,flip_signature  returns -1 if $x_i<x_j$ and O else. 
It will return -2 in case the Graph G has a loop. 
"""
function flip_signature(G::FeynmanGraph ,p::Vector{Int64},a::Vector{Int64}) #graph G, list p, branch type a
    ee=edges(G)
    b=zeros(Int,length(a))
    for (i, (ai, ev)) in enumerate(zip(a, ee))
        if ai == 0 && src(ev) != dst(ev)
            if preimg(p, src(ev)) < preimg(p, dst(ev))
                b[i] = -1
            else
                b[i] = 0
            end
        elseif  src(ev) == dst(ev)
            b[i] = -2
        elseif ai != 0 && src(ev) != dst(ev)
            b[i] = ai
        end
    end
    return b
end
function signature_and_multiplicities_order(G::FeynmanGraph, a::Vector{Int64},o::Vector)
    b=Vector{Tuple{Int64, Vector{Int64}}}()
    push!(b,(1,flip_signature(G,o,a))) 
return b
end
# flip_signature regroup all Orders with the same signature so same
@doc raw"""
    signature_and_multiplicities( G::FeynmanGraph, a::Vector{Int64})

 returns flip_signature and their multiplicities.
# Examples
```julia 
julia> using GromovWitten

julia> G=FeynmanGraph([(1, 1), (1, 2), (2, 3), (3, 1)])
FeynmanGraph([(1, 1), (1, 2), (2, 3), (3, 1)])

julia> a=[2,0,0,1];

julia> signature_and_multiplicities(G,a)
2-element Vector{Tuple{Int64, Vector{Int64}}}:
 (4, [-2, -1, 0, 1])
 (2, [-2, 0, 0, 1])
```
"""
function signature_and_multiplicities(G::FeynmanGraph, a::Vector{Int64})
    ee = edges(G)
    p = Vector{Int64}()
    b = Vector{Tuple{Int64, Vector{Int64}}}()
    l = zeros(Int, nv(G))
    y = Vector{Vector{Int64}}()
    if count_zero(a)<=1
        push!(b,(factorial( nv(G)) ,a))
         return b
    else
         for (ev,ai) in zip(ee,a)
            if ai==0 && src(ev) != dst(ev)
                l[src(ev)] =1
                l[dst(ev)] =1
            end
        end
                #println(l)

        for (i,li) in enumerate(l)
            if li==1
                push!(p,i)
            end
        end 
        p=collect(permutations(p))

        for ga in p 
            fl=flip_signature(G,ga,a)
            push!(y,fl)
        end
        dd=div(factorial( nv(G) ) , length( p ) )
        py=countmap(y)
        for (key, val) in py
            push!(b, (dd*val,key ))
        end
        if length(b) == 1
            return b
        else
            group = Vector{Tuple{Int64, Vector{Int64}}}()

            for (n, values1) in b
                mm = 2 * n
                if (n, values1) in group || (mm, values1) in group
                    continue
                end

                equiv = false

                for (m, values2) in b
                    if (m, values2) in group || (2 * m, values2) in group
                        continue
                    end

                    if n == m && values2 == replace(values1)
                        equiv = true
                        break
                    end
                end

                mn = 2 * n

                if equiv
                    push!(group, (mm, values1))
                end
            end

            # Convert Set to Vector for consistency with the original return type
            return group
        end
    end
end


function find_equal_pairs(ve::Vector{Edge})
    equal_pairs = Dict{Edge, Vector{Int}}()
    for (i, pair) in enumerate(ve)
        if haskey(equal_pairs, pair)
            push!(equal_pairs[pair], i)
        else
            equal_pairs[pair] = [i]
        end
    end
    indices = [v for v in values(equal_pairs) if length(v) > 1]
    return indices
end
function vector_to_monomial(F::FeynmanIntegral,v::Vector{Int64})
    S=polynomial_ring(QQ, :x =>1:nv(F.G), :q =>1:ne(F.G), :z => 1:nv(F.G))
    q=S[3]
    v = 2 * v
    vec=[]
    for i in eachindex(v)
        push!(vec,q[i]^v[i])
    end
    poly = prod(vec)
    return poly
end
function generate_permutation(l::Vector{Int64}, indices::Vector{Vector{Int64}})
    original_l = copy(l)
    permuted_lists = Set{Vector{Int64}}()  # Use a Set to ensure uniqueness
    group_permutations = []

    for group_indices in indices
        group_elements = l[group_indices]
        push!(group_permutations, collect(permutations(group_elements)))
    end

    # Generate permutations for each group
    for permuted_indices in Iterators.product(group_permutations...)
        temp_l = copy(original_l)
        for (i, group_indices) in enumerate(indices)
            for (j, idx) in enumerate(group_indices)
                temp_l[idx] = permuted_indices[i][j]
            end
        end
        push!(permuted_lists, temp_l)
    end

    return collect(permuted_lists)
end