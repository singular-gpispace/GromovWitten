module tropicalfeynman

# Write your package code here.
using Oscar, Graphs, Combinatorics,StatsBase

#ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
 struct graphe
    edge::Vector
end
import Graphs: nv,ne,Edge,dst, src


export graphe,nv,ne,PolynomialRing,QQ
export polynomialring,constterm, proterm, propagator, coefterm, partition, preimg, sgn,flip,
specificFeynmanIntegral, feynmanIntegral,feynmanIntegralSum,subt,polynomialringV,consttermV, protermV, coefterm2Z, coeftermQ,coeftermX, sgnV,flipV,
specificFeynmanIntegralV, feynmanIntegralV,feynmanIntegralSumV,subtV


include("graphe.jl")
include("coeftermV.jl")
include("coefterm.jl")
include("ceoftermGV.jl")
include("feynmanIntegral.jl")
include("propagator.jl")

end
