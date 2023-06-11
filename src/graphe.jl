function polynomialring(G::graphe, x::String, q::String, z::String)
    return @polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)], z[1:nv(G)])
end
function polynomialring(G::graphe, x::String, q::String)
    return @polynomial_ring(QQ, x[1:nv(G)], q[1:ne(G)])
end
edge(G::graphe)=Edge.(G.edge)
nv(G::graphe)=nv(DiGraph(Edge.(G.edge)))
ne(G::graphe)=length(Edge.(G.edge))