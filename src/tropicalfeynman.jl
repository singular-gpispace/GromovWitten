module tropicalfeynman

# Write your package code here.
using Oscar, Graphs, Combinatorics,StatsBase

el = Edge.([ (2, 4), (4,1),(4, 1), (1, 3) ,(2, 3),(2,3)])
ee=collect(el)
G= DiGraph(el)

R, x, q = PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:ne(G))
export 
constterm, term, propagator,coefterm,partition, preimg, sgn,flip,
specificFeynmanIntegral, feynmanIntegral,feynmarenIntegralParallel,feynmanIntegralSum,sub



include("coefterm.jl")
include("feynmanIntegral.jl")
#include("propagator.jl")

end
