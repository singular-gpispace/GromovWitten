struct graphe
    edge::Vector
end
function polynomialring(G::graphe)
    return  R, x, q = PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:ne(G))
 
end
function polynomialringV(G::graphe)
    return  PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:ne(G),"z" => 1:nv(G))
end
edge(G::graphe)=Edge.(G.edge)
nv(G::graphe)=nv(DiGraph(Edge.(G.edge)))
ne(G::graphe)=length(Edge.(G.edge))
#R, x, q = PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:ne(G))

