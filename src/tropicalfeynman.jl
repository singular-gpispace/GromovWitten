@doc raw"""
 tropicalfeynman  is a package for computing Feynman Integral via Gromov-Witten invariant.
"""
module tropicalfeynman
# Write your package code here.
using Nemo, Graphs, Combinatorics,StatsBase

import Graphs: nv,ne,Edge,dst, src
export graph,nv,ne,PolynomialRing,QQ
export polynomialring,constterm, proterm, propagator, coefterm, partition, preimg, sgn,flip,flipo,
feynman_integral_branchtype, feynman_integral_degree,feynman_integral_degreeSum,subt,consttermV, protermV, coefterm2Z, coeftermQ,coeftermX, sgnV,flipV,
feynman_integral_branchtype_order,feynman_integral_degree_order,lis,filter_term

include("graph.jl")
include("coeftermV.jl")
include("coefterm.jl")
include("feynmanIntegral.jl")
include("propagator.jl")
end
