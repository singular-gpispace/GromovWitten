###############################################################################
#                                                                             #
#            graph.jl : define graph and Polynomial Ring.                     #
#                                                                             #
###############################################################################
struct FeynmanGraph
    edge::Vector{Tuple{Int, Int}}
end

nv(G::FeynmanGraph) = length(union([e[1] for e in G.edge], [e[2] for e in G.edge]))
ne(G::FeynmanGraph) = length(G.edge)

function feynman_graph(edges::Vector{Tuple{Int, Int}})
    return FeynmanGraph(edges)
end
#=function polynomial_ring(G::FeynmanGraph, x::String, q::String, z::String)
    ee=Edge.(G.edge)
    if length(ee)==0
        throw(DomainError(G,"G must be non empty graph"))
    else
        return @polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)], z[1:nv(G)])
    end
end
function polynomial_ring(G::FeynmanGraph, x::String, q::String)
    ee=Edge.(G.edge)
    if length(ee)==0
        throw(DomainError(G,"G must be non empty graph"))
    else
        return polynomial_ring(QQ, :x =>1:nv(G), :q =>1:ne(G), :z => 1:nv(G))
    end
end=#
const FeynmanRing{T} =Tuple{Ring,Vararg{Vector{T}}} where T <:QQMPolyRingElem
 
struct FeynmanIntegral
    G::FeynmanGraph
    integral_cache::Dict{Symbol, Dict{Vector{Int64}, QQMPolyRingElem}}
    S::FeynmanRing

    #External constructor if S is given 
    function FeynmanIntegral(G::FeynmanGraph, S::FeynmanRing)
        return new(G, Dict{Symbol, Dict{Vector{Int64}, QQMPolyRingElem}}(), S)
    end
    
    # Inner constructor with default values for S
    function FeynmanIntegral(G::FeynmanGraph)
        S=polynomial_ring(QQ, :x =>1:nv(G), :q =>1:ne(G), :z => 1:nv(G))
        return new(G, Dict{Symbol, Dict{Vector{Int64}, QQMPolyRingElem}}(), S)
    end
    function FeynmanIntegral(ve::Vector{Tuple{Int, Int}})
        G = FeynmanGraph(ve)
        F = FeynmanIntegral(G)
        return  F
    end
end
function feynman_integral(G::FeynmanGraph)
    return FeynmanIntegral(G)
end
 edge(F::FeynmanIntegral) = Edge.(F.G.edge)