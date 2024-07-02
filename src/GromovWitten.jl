@doc raw"""
GromovWitten  is a package for computing Gromov-Witten invariant via Feynman Integral  .
"""
module GromovWitten

using Nemo

import Combinatorics: permutations
export next_partition, combination
export count_member
export cache_integral_result
export coefterm
export constterm
export eisenstein_series
export express_as_eisenstein_series
export express_as_powers
export FeynmanGraph
export FeynmanIntegral
export feynman_graph
export feynman_integral
export feynman_integral_branch_type
export feynman_integral_branch_type_cache
export feynman_integral_branch_type_order
export feynman_integral_deg
export feynman_integral_degree
export feynman_integral_degree_cache
export feynman_integral_degree_order
export feynman_integral_degree_sum
export feynman_integral_degree_sum_cache
export feynman_integral_degree_sum_order
export filter_term
export filter_vector
export find_equal_pairs
export flip
export flip_signature
export generate_permutation
export get_integral_from_cache
export inv_sfunction
export lis
export loopterm
export matrix_of_integral
export ne, src
export nv, dst
export edges, edg
export partition
export polynomial_ring
export polynomial_to_matrix
export preimg
export proterm
export QQ
export quasimodular_form, quasimodularity_form
export quasi_matrix, quasimodular_matrix
export replace
export sfunction
export signature_and_multiplicities
export signature_and_multiplicities_order
export solve_polynomial_system
export substitute
export sum_of_coeff
export sum_of_divisor_powers
export vector_to_monomial
export number_of_monomials
export compos_iterate
include("graph.jl")
include("coeftermV.jl")
include("coefterm.jl")
include("feynmanIntegral.jl")
include("propagator.jl")
include("quasimodular.jl")
include("hashtable.jl")

end
