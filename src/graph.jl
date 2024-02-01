###############################################################################
#                                                                             #
#            graph.jl : define graph and Polynomial Ring.                     #
#                                                                             #
###############################################################################
struct Edge
    src::Int
    dst::Int
end
Edge(e::Tuple{Int, Int}) = Edge(e[1], e[2])
src(e::Edge) = e.src
dst(e::Edge) = e.dst

Base.show(io::IO, e::Edge) = show(io, (src(e), dst(e)))

struct FeynmanGraph
    _edges::Vector{Edge}
end
FeynmanGraph(edges::Vector{Tuple{Int, Int}}) = FeynmanGraph(Edge.(edges))

edges(G::FeynmanGraph) = G._edges
nv(G::FeynmanGraph) = length(union([e.src for e in G._edges], [e.dst for e in G._edges]))
ne(G::FeynmanGraph) = length(G._edges)

Base.show(io::IO, G::FeynmanGraph) = print(io, "FeynmanGraph(", [ (src(e), dst(e)) for e in edges(G) ], ")")

function feynman_graph(edges::Vector{Tuple{Int, Int}})
    return FeynmanGraph(edges)
end
#=function polynomial_ring(G::FeynmanGraph, x::String, q::String, z::String)
    ee=edges(G)
    if length(ee)==0
        throw(DomainError(G,"G must be non empty graph"))
    else
        return @polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)], z[1:nv(G)])
    end
end
function polynomial_ring(G::FeynmanGraph, x::String, q::String)
    ee=edges(G)
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
