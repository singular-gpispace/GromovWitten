@doc raw"""
 tropicalfeynman  is a package for computing Feynman Integral via Gromov-Witten invariant.

 For more information see
"""
module tropicalfeynman
# Write your package code here.
using Nemo, Graphs, Combinatorics,StatsBase

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
