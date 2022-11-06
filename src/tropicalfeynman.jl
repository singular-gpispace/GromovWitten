module tropicalfeynman

# Write your package code here.
using Oscar, Graphs, Combinatorics,StatsBase

el = Edge.([(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)])
ee=collect(el)
G= DiGraph(el)

R, x, q = PolynomialRing(QQ, "x" => 1:nv(G), "q" => 1:length(ee))
a=[1,2,3,2,0,1]
export 
constterm, term, propagator,coefterm,partition, preimg, sgn,flip,
specificFeynmanIntegral, feynmanIntegral,feynmarenIntegralParallel,feynmanIntegralSum,sub


include("coefterm.jl")
include("feynmanIntegral.jl")
include("propagator.jl")

end
