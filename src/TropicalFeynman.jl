@doc raw"""
 tropicalfeynman  is a package for computing Feynman Integral via Gromov-Witten invariant.
"""
module TropicalFeynman
# Write your package code here.
using Nemo, Graphs, Combinatorics,StatsBase

import Graphs: nv,ne,Edge,dst, src
export graph,nv,ne,PolynomialRing,QQ
export polynomial_ring,constterm, proterm, propagator, coefterm, partition, preimg, substitute,flip,signature_and_multiplicities_order,
feynman_integral_branchtype, feynman_integral_degree,feynman_integral_degree_sum,subt,consttermV, protermV,  flip_signature,signature_and_multiplicities,
feynman_integral_branchtype_order,feynman_integral_degree_order,lis,filter_term,filter_vector,
sum_of_divisor_powers,express_as_eisenstein_series,express_as_powers,polynomial_to_matrix,matrix_of_integral,
solve_polynomial_system,quasi_matrix

include("graph.jl")
include("coeftermV.jl")
include("coefterm.jl")
include("feynmanIntegral.jl")
include("propagator.jl")
include("quasimodular.jl")
end
