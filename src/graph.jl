###############################################################################
#                                                                             #
#           graph.jl : define graph and Polynomial Ring.                     #
#                                                                             #
###############################################################################

struct graph
    edge::Vector
end
edge(G::graph)=Edge.(G.edge)
nv(G::graph)=nv(DiGraph(Edge.(G.edge)))
ne(G::graph)=length(Edge.(G.edge))
function polynomial_ring(G::graph, x::String, q::String, z::String)
    return @polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)], z[1:nv(G)])
end
function polynomial_ring(G::graph, x::String, q::String)
    return @polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)])
end
