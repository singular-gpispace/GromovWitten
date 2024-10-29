julia --project


using GromovWitten
ve = [(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)]
ve = [(1, 2), (1, 2), (1, 3), (2, 4), (3, 4), (3, 5), (4, 6), (5, 6), (5, 6)]
ve = [(1, 2), (1, 2), (1, 2)]

#ve = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
#ve = [(1, 2), (1, 2), (1, 3), (2, 4), (3, 4), (3, 5), (4, 6), (5, 6), (5, 7), (6, 8), (7, 8), (7, 8)];
#ve = [(1, 2), (1, 3), (1, 2), (2, 4), (3, 4), (3, 4)]
#ve = [(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 5), (4, 6), (5, 6), (5, 6)]

F = FeynmanIntegral(ve)
weightmax=6
m=number_of_monomials(6)
 Iq=substitute(feynman_integral_degree_sum(F, m))
 quasimodular_form(Iq,weightmax)