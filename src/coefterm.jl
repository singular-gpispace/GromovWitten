###############################################################################
#                                                                             #
#     coefterm.jl : compute differents  coefficients of polynomial.           #
#                                                                             #
###############################################################################

#  coefterm compute the coefficients term x1^0,...,xn^0, By I(q)=coeff[x1^0,...,xn^0](P(x,q))
#Where I(q) is the Feynman Integral and P(x,q) the Propagator.
function coefterm(x::Vector,q::Vector,G::graph ,p::fmpq_mpoly,d::Integer;l=zeros(Int,nv(G)))
   #here l is leak vector of the graph G.
    ee=Edge.(G.edge) 
    G=DiGraph(Edge.(G.edge)) # convert from graph to Graph (so we can use nv(G))
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
end 
@doc raw"""
partition(k::Integer, n::Integer)    

**Note**:This function returns the number of partitions of $n$ into fixed  $k$ parts. 
# Examples
```jldoctest
julia> partition(3,4)
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
````
"""
function partition(n::Integer, k::Integer)
    if k == 1
        return [[n]]
    end
    if n == 0
        return [[0] * k]
    end
    result = Array{Int64}[]
    for A in with_replacement_combinations(1:n, k)
        item = zeros(Int64, n)
        for a in A
            item[a] += 1
        end
        push!(result, item)
    end
    return result
end
#give the position of the vertices xi in the list L. 
function preimg(L::Vector, xi::Integer)
    for (i,Li) in enumerate(L)
        if Li==xi
            return i # i is then the position of xi in L
            break
        end
    end
end

# This first signature determine the flip signature. Ω=[x1,...,xn] is a given Order 
#and a is a branche type. It returns -1 if xi<xj and O else. 
function sgn(G::graph ,Ω::Vector,a::Vector) #graph G, list p, branch type a  
    ee=Edge.(G.edge)
    b=zeros(Int,length(a))
    for (i,(ai,ev)) in enumerate(zip(a,ee))
       if ai==0
           if preimg(Ω,src(ev))<preimg(Ω,dst(ev))
                    b[i]=-1
                else 
                    b[i]=0
            end
            else
                b[i]=ai # So only O entry of a are affected. 
        end         
    end
            return b
end
# Same like before, here it returns -2 in case the Graph G has a loop. 
# It detects the loop in the graph. 
function sgnV(G::graph ,p::Vector,a::Vector) #graph G, list p, branch type a
    ee=Edge.(G.edge)
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

function flip( G::graph,  a::Vector)

    ee=Edge.(G.edge)
    p=[]
    b=Vector{Vector{Any}}()

    l=zeros(Int,nv(G))
    y=Vector{Any}()

    for (ev,ai) in zip(ee,a)
        if ai==0
            l[src(ev)] =1
            l[dst(ev)] =1
        end
    end
    for (i,li) in enumerate(l)
           if li==1
            push!(p,i)
        end
    end   
    p=collect(permutations(p))

    for ga in p 
        push!(y,sgn(G,ga,a))    
    end
    dd=div(factorial( nv(G) ) , length( p ) )
 
    for (key, val) in countmap(y)
        push!(b,[ val*dd, key]) 
    end
    return b   
end
# Flip signature regroup all Orders with the same signature so same
function flipV( G::graph, a::Vector)

    ee=Edge.(G.edge)

    p=[]
    b=Vector{Vector{Any}}()

    l=zeros(Int,nv(G))
    y=Vector{Any}()

    for (ev,ai) in zip(ee,a)
        if ai==0 && src(ev) != dst(ev)
            l[src(ev)] =1
            l[dst(ev)] =1
        end
    end

    for (i,li) in enumerate(l)
        if li==1
            push!(p,i)
    end

    end 
    p=collect(permutations(p))

    for ga in p 
        push!(y,sgnV(G,ga,a)) 
    end
    dd=div(factorial( nv(G) ) , length( p ) )

    for (key, val) in countmap(y)
        push!(b,[ val*dd, key]) 
    end
    return b 
end
function flipo(G::graph, a::Vector,o::Vector)
    b=Vector{Vector{Any}}()
    y=Vector{Any}()
    push!(y,sgnV(G,o,a)) 
    for (key, val) in countmap(y)
        push!(b,[ val, key]) 
    end
    return b
end