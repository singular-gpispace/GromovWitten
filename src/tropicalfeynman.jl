module tropicalfeynman
# Write your package code here.
using Nemo, Graphs, Combinatorics,StatsBase

#ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]
import Graphs: nv,ne,Edge,dst, src
export graphe,nv,ne,PolynomialRing,QQ
export polynomialring,constterm, proterm, propagator, coefterm, partition, preimg, sgn,flip,flipo,
specificFeynmanIntegral, feynmanIntegral,feynmanIntegralSum,subt,consttermV, protermV, coefterm2Z, coeftermQ,coeftermX, sgnV,flipV,
specificFeynmanIntegralo,feynmanIntegralo

include("graphe.jl")
include("coeftermV.jl")
include("coefterm.jl")
include("feynmanIntegral.jl")
include("propagator.jl")
end
