@doc raw"""
GromovWitten  is a package for computing Gromov-Witten invariant via Feynman Integral  .
"""
module GromovWitten
# Write your package code here.
using Nemo, Graphs, Combinatorics,StatsBase
export  polynomial_ring
import Graphs: nv,ne,Edge,dst, src
export FeynmanGraph,nv,ne,PolynomialRing,QQ,FeynmanIntegral
export polynomial_ring,constterm, proterm, propagator, coefterm, partition, preimg, substitute,flip,signature_and_multiplicities_order,
feynman_integral_branch_type, feynman_integral_degree,feynman_integral_degree_sum_order,feynman_integral_degree_sum,sum_of_coeff,loopterm,  flip_signature,signature_and_multiplicities,
feynman_integral_branch_type_order,feynman_integral_degree_order,lis,filter_term,filter_vector,
sum_of_divisor_powers,express_as_eisenstein_series,express_as_powers,polynomial_to_matrix,matrix_of_integral,
solve_polynomial_system, quasi_matrix,sfunction,inv_sfunction
export eisenstein_series,quasimodular_form
export cache_integral_result, get_integral_from_cache, feynman_integral_branch_type_cache, feynman_integral_degree_cache,feynman_integral_degree_sum_cache
export feynman_integral,generate_permutation,find_equal_pairs,vector_to_monomial,replace,feynman_integral_deg
include("graph.jl")
include("coeftermV.jl")
include("coefterm.jl")
include("feynmanIntegral.jl")
include("propagator.jl")
#include("quasimodular.jl")
include("hashtable.jl")
end
