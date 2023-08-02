###############################################################################
#                                                                             #
#            graph.jl : define graph and Polynomial Ring.                     #
#                                                                             #
###############################################################################
struct FeynmanGraph
    edge::Vector
end
 # edge(G::FeynmanGraph) = Edge.(G.edge)
    nv(G::FeynmanGraph) = nv(DiGraph(Edge.(G.edge)))
    ne(G::FeynmanGraph) = length(Edge.(G.edge))

function polynomial_ring(G::FeynmanGraph, x::String, q::String, z::String)
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
        return @polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)])
    end
end
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
        S=@polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)], z[1:nv(G)])
        return new(G, Dict{Symbol, Dict{Vector{Int64}, QQMPolyRingElem}}(), S)
    end
end